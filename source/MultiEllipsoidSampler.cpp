#include "MultiEllipsoidSampler.h"

// MultiEllipsoidSampler::MultiEllipsoidSampler()
//
// PURPOSE: 
//      Class constructor.
//
// INPUT:
//      printOnTheScreen:               true if the results are to be printed on the screen, false otherwise
//      ptrPriors:                      vector of Prior class objects containing the priors used in the problem.
//      likelihood:                     Likelihood class object used for likelihood sampling.
//      metric:                         Metric class object to contain the metric used in the problem.
//      clusterer:                      Clusterer class object specifying the type of clustering algorithm to be used.
//      initialNlivePoints:                Initial number of active points to start the nesting process
//      NlivePointsMinimum:                Minimum number of active points allowed in the nesting process
//      initialEnlargementFraction:     Initial value of the enlargement for the ellipsoids
//      shrinkingRate:                  Shrinking rate of the enlargement factors, between 0 and 1.
//

MultiEllipsoidSampler::MultiEllipsoidSampler(const bool printOnTheScreen, vector<Prior*> ptrPriors, 
                                             Likelihood &likelihood, Metric &metric, Clusterer &clusterer,
                                             const int initialNlivePoints, const int minNlivePoints, 
                                             const double initialEnlargementFraction, const double shrinkingRate)
: NestedSampler(printOnTheScreen, initialNlivePoints, minNlivePoints, ptrPriors, likelihood, metric, clusterer),
  ellipsoidMatrixDecompositionIsSuccessful(true),
  initialEnlargementFraction(initialEnlargementFraction),
  shrinkingRate(shrinkingRate),
  uniform(0.0, 1.0)
{
}










// MultiEllipsoidSampler::~MultiEllipsoidSampler()
//
// PURPOSE: 
//      Base class destructor.
//

MultiEllipsoidSampler::~MultiEllipsoidSampler()
{

}











// MultiEllipsoidSampler::drawWithConstraint()
//
// PURPOSE:
//      Draws a new point from one of the ellipsoids built from the sample of identifed clusters.
//      The drawing process follows the adopted prior distributions. If no new point satisfying the 
//      hard likelihood constraint is found, then a false boolean value is returned. 
//      The ellipsoid to draw from is randomly selected according to its volume (Feroz F., Hobson M. P., 2008, MNRAS, 384, 449). 
//      In case ellipsoids are overlapping, if the drawn point falls in a region where ellipsoids overlap, 
//      then the point is selected with a probability inverse to the number of ellipsoids overlapping in that region.
//
// INPUT:
//      totalSample:                Eigen Array matrix of size (Ndimensions, NlivePoints)
//                                  containing the total sample of active points at a given nesting iteration
//      Nclusters:                  Optimal number of clusters found by clustering algorithm
//      clusterIndices:             Indices of clusters for each point of the sample
//      clusterSizes:               A vector of integers containing the number of points belonging to each cluster
//      drawnPoint:                 Eigen Array matrix of size (Ndimensions,Ndraws) to contain the
//                                  coordinates of the drawn point to be used for the next nesting loop. 
//                                  When used for the first time, the array contains the coordinates of the 
//                                  worst nested object to be updated.
//      logLikelihoodOfDrawnPoint:  the log(likelihood) value of the new drawn point, to be stored in the global
//                                  array containing the log(likelihood) values for each nested iteration.
//      maxNdrawAttempts:           Maximum number of attempts allowed when drawing from a single ellipsoid.
//      
//
// OUTPUT:
//      A boolean value that is true if a new point in the sampling process is found and false otherwise.
//

bool MultiEllipsoidSampler::drawWithConstraint(const RefArrayXXd totalSample, const unsigned int Nclusters, const vector<int> &clusterIndices,
                                               const vector<int> &clusterSizes, RefArrayXd drawnPoint, 
                                               double &logLikelihoodOfDrawnPoint, const int maxNdrawAttempts)
{    
    assert(totalSample.cols() == clusterIndices.size());
    assert(drawnPoint.size() == totalSample.rows());
    assert(Nclusters > 0);


    // Compute the ellipsoids corresponding to the clusters found by the clustering algorithm.
    // This involves computing the barycenter, covariance matrix, eigenvalues and eigenvectors
    // for each ellipsoid/cluster.

    computeEllipsoids(totalSample, Nclusters, clusterIndices, clusterSizes);


    // Find which ellipsoids are overlapping and which are not
    
    vector<unordered_set<int>> overlappingEllipsoidsIndices;
    findOverlappingEllipsoids(overlappingEllipsoidsIndices);


    // If the ellipsoid matrix decomposition fails, return to main nested sampling loop with no new point drawn

    if (!ellipsoidMatrixDecompositionIsSuccessful)
    {
        return false;
    }


    // Get the hyper-volume for each of the ellipsoids and normalize it 
    // to the sum of the hyper-volumes over all the ellipsoids

    vector<double> normalizedHyperVolumes(Nellipsoids);
    
    for (int n=0; n < Nellipsoids; ++n)
    {
        normalizedHyperVolumes[n] = ellipsoids[n].getHyperVolume();
    }

    double sumOfHyperVolumes = accumulate(normalizedHyperVolumes.begin(), normalizedHyperVolumes.end(), 0.0, plus<double>());

    for (int n=0; n < Nellipsoids; ++n)
    {
        normalizedHyperVolumes[n] /= sumOfHyperVolumes;
    }


    // Select an ellipsoid. Sample a point from the prior within this ellipsoid, such that its Likelihood  
    // value is better than the worst likelihood of all points currently contained in this ellipsoid.

    // Pick an ellipsoid with a probability according to its normalized hyper-volume
    // First generate a uniform random number between 0 and 1

    double uniformNumber = uniform(engine);


    // Select the ellipsoid that makes the cumulative hyper-volume greater than this random
    // number. Those ellipsoids with a larger hyper-volume will have a greater probability to 
    // be chosen.

    double cumulativeHyperVolume = normalizedHyperVolumes[0];
    int indexOfSelectedEllipsoid = 0;
    
    while (cumulativeHyperVolume < uniformNumber)
    {
        indexOfSelectedEllipsoid++;
        cumulativeHyperVolume += normalizedHyperVolumes[indexOfSelectedEllipsoid];
    }


    // Draw a new point in this Ellipsoid, but with the constraints mentioned above.

    bool newPointIsFound = false;

    
    // Drawing a new point with the constraints mentioned above may prove difficult.
    // Hence, we won't try more than 'maxNdrawAttempts'.

    int NdrawAttempts = 0;

    while ((newPointIsFound == false) & (NdrawAttempts < maxNdrawAttempts))
    {
        // Keep count of the number of attempts

        NdrawAttempts++;


        // Draw a new point inside the ellipsoid

        ellipsoids[indexOfSelectedEllipsoid].drawPoint(drawnPoint);


        // Check if the new point is also in other ellipsoids. If the point happens to be 
        // in N overlapping ellipsoids, then accept it only with a probability 1/N. If we
        // wouldn't do this, the overlapping regions in the ellipsoids would be oversampled.

        if (!overlappingEllipsoidsIndices[indexOfSelectedEllipsoid].empty())
        {
            // There are overlaps, so count the number of ellipsoids to which the new
            // point belongs
            
            int NenclosingEllipsoids = 1;

            for (auto index = overlappingEllipsoidsIndices[indexOfSelectedEllipsoid].begin();
                      index != overlappingEllipsoidsIndices[indexOfSelectedEllipsoid].end();
                      ++index)
            {
                if (ellipsoids[*index].containsPoint(drawnPoint))  NenclosingEllipsoids++;
            }


            // Only accept the new point with a probability = 1/NenclosingEllipsoids. 
            // If it's not accepted, go immediately back to the beginning of the while loop, 
            // and draw a new point inside the ellipsoid.

            uniformNumber = uniform(engine);
            newPointIsFound = (uniformNumber < 1./NenclosingEllipsoids);
            if (!newPointIsFound) continue;
        }
        else
        {
            // There are no ellipsoids overlapping with the selected one, so the point
            // is automatically accepted

            newPointIsFound = true;
        }


        // The point should not only be drawn inside the ellipsoid, it should also be drawn
        // from the prior. Therefore, accept the point only with the probability given by the
        // prior, so that the regions inside the ellipsoid with a higher prior density will 
        // be sampled more than the regions with a lower prior density. 
        // Since different coordinates of our new point may have different priors, 
        // we need to check this for all the priors.
        
        int beginIndex = 0;

        for (int priorIndex = 0; priorIndex < ptrPriors.size(); ++priorIndex)
        {
            // Figure out the number of parameters (=coordinates) that the current prior covers

            const int NdimensionsOfPrior = ptrPriors[priorIndex]->getNdimensions();


            // Define a subset of the new point, consisting of those coordinates covered by the 
            // same (current) prior distribution.

            ArrayXd subsetOfNewPoint = drawnPoint.segment(beginIndex, NdimensionsOfPrior);


            // Check if the new point is accepted according to the corresponding prior distribution.
            // If not, exit the loop and draw a new point from the beginning.

            newPointIsFound = ptrPriors[priorIndex]->drawnPointIsAccepted(subsetOfNewPoint);
                
            if (!newPointIsFound)
            break;


            // Move the beginIndex on to the next set of coordinates covered by the prior

            beginIndex += NdimensionsOfPrior;
        }

        if (!newPointIsFound)
        {
            // The new point failed the prior criterion, so go back to the start of the while loop
            // and draw a new point inside the selected ellipsoid
        
            continue;
        }


        // Finally, the point should not only be drawn inside the ellipsoid and according to the prior
        // density, but it should also have a likelihood that is larger than the one of the worst point.
        // We check this criterion only after the prior criterion, because often the likelihood is
        // much more time consuming to compute than the prior.

        logLikelihoodOfDrawnPoint = likelihood.logValue(drawnPoint);

        if (logLikelihoodOfDrawnPoint < worstLiveLogLikelihood)
        {
            // The new point does not fulfill the likelihood criterion. Flag it as such,
            // and go back to the start of the while loop, and draw a new point inside the
            // ellipsoid.

            newPointIsFound = false;
        }

    } // end while-loop (newPointIsFound == false)


    // Depending on whether we found a new point or not, return true or false.

    if (newPointIsFound)
    {
        return true;
    }
    else
    {
        return false;
    }
}











// MultiEllipsoidSampler::verifySamplerStatus()
//
// PURPOSE:
//      Verifies whether the status of the sampler in use is successful.
//      If not, it prints the related error message and returns a false boolean value.
//      In this case, it concerns the error caused by a failure in the ellipsoid matrix 
//      decomposition.
//
// OUTPUT:
//      void
//

bool MultiEllipsoidSampler::verifySamplerStatus()
{
    if (!ellipsoidMatrixDecompositionIsSuccessful)
    {
        cout << "Ellipsoid Matrix decomposition failed." << endl;
        cout << "Quitting program." << endl;
        return false;
    }
    else
    {
        return true;
    }
}











// MultiEllipsoidSampler::computeEllipsoids()
//
// PURPOSE:
//      Computes covariance matrices and center coordinates of all the ellipsoids
//      associated to each cluster of the sample. The eigenvalue decomposition
//      is also done and eigenvalues are enlarged afterwards. All the results
//      are stored in the private data members.
//
// INPUT:
//      totalSample(Ndimensions, NlivePoints):     Complete sample (spread over all clusters) of points
//      Nclusters:                              The number of clusters identified by the clustering algorithm
//      clusterIndices(NlivePoints):               For each point, the integer index of the cluster to which it belongs
//      clusterSizes(Nclusters):                A vector of integers containing the number of points belonging to each cluster
//
// OUTPUT:
//      void
//

void MultiEllipsoidSampler::computeEllipsoids(RefArrayXXd const totalSample, const unsigned int Nclusters, 
                                              const vector<int> &clusterIndices, const vector<int> &clusterSizes)
{
    assert(totalSample.cols() == clusterIndices.size());
    assert(totalSample.cols() >= Ndimensions + 1);            // At least Ndimensions + 1 points are required to start.


    // Compute "sorted indices" such that clusterIndices[sortedindices[k]] <= clusterIndices[sortedIndices[k+1]]

    vector<int> sortedIndices = Functions::argsort(clusterIndices);


    // beginIndex will take values such that the indices for one particular cluster (# n) will be in 
    // sortedIndex[beginIndex, ..., beginIndex + clusterSize[n] - 1]      

    int beginIndex = 0;


    // Clear whatever was in the ellipsoids collection

    ellipsoids.clear();


    // Create an Ellipsoid for each cluster (provided it's large enough)

    for (int i = 0; i < Nclusters; i++)
    {   
        // Skip cluster if number of points is not large enough

        if (clusterSizes[i] < Ndimensions + 1) 
        {
            // Move the beginIndex up to the next cluster

            beginIndex += clusterSizes[i];


            // Continue with the next cluster

            continue;
        }
        else
        {
            // The cluster is indeed large enough to compute an Ellipsoid.

            // Copy those points that belong to the current cluster in a separate Array
            // This is because Ellipsoid needs a contiguous array of points.

            ArrayXXd sampleOfOneCluster(Ndimensions, clusterSizes[i]);

            for (int n = 0; n < clusterSizes[i]; ++n)
            {
                sampleOfOneCluster.col(n) = totalSample.col(sortedIndices[beginIndex+n]);
            }
    

            // Move the beginIndex up to the next cluster

            beginIndex += clusterSizes[i];


            // Compute the new enlargement fraction (it is the fraction by which each axis of an ellipsoid is enlarged).
            // This allows for improving the efficiency of the sampling by increasing the chance of having more
            // points of the cluster falling inside the bounding ellipsoid.

            double enlargementFraction = updateEnlargementFraction(clusterSizes[i]);
           

            // Add ellipsoid at the end of our vector

            ellipsoids.push_back(Ellipsoid(sampleOfOneCluster, enlargementFraction));
        }
    }

    Nellipsoids = ellipsoids.size();
}











// MultiEllipsoidSampler::findOverlappingEllipsoids()
//
// PURPOSE:
//      For each ellipsoid, determine which other ellipsoids overlap with it, and
//      store the indices of those overlapping ellipsoids.
//
// INPUT:
//      overlappingEllipsoidsIndices[0..Nellipsoids-1]:     a vector of unordered_set of integers
//                                                          whose elements will contain the indices
//                                                          corresponding to the overlapping ellipsoids
//
// OUTPUT:
//      overlappingEllipsoidsIndices will have changed. Element [i] contains an unordered_set<> with 
//      the indices of all ellipsoids that overlap with ellipsoid #i
//

void MultiEllipsoidSampler::findOverlappingEllipsoids(vector<unordered_set<int>> &overlappingEllipsoidsIndices)
{
    // Remove whatever was in the container before

    overlappingEllipsoidsIndices.clear();
   

    // Make sure that the indices container has the right size

    overlappingEllipsoidsIndices.resize(ellipsoids.size());


    // If Ellipsoid i overlaps with ellipsoid j, then ellipsoid j also overlaps with i.
    // The indices are kept in an unordered_set<> which automatically gets rid of duplicates.
    // If an eigenvalues decomposition error occurs, it is stored in the 
    // boolean variable ellipsoidMatrixDecompositionIsSuccessful

    for (int i = 0; i < Nellipsoids-1; ++i)
    {
        for (int j = i+1; j < Nellipsoids; ++j)
        {
            if (ellipsoids[i].overlapsWith(ellipsoids[j], ellipsoidMatrixDecompositionIsSuccessful))
            {
                overlappingEllipsoidsIndices[i].insert(j);
                overlappingEllipsoidsIndices[j].insert(i);
            }
        }
    }
}











// MultiEllipsoidSampler::updateEnlargementFraction()
//
// PURPOSE:
//      Updates the enlargementFraction adopted for the axes of the ellipsoid.
//      The formula takes into account the number of dimensions of the problem
//      and it is a modified version of the one adopted by Feroz F. et al. 2008.
//
// INPUT:
//      clusterSize:    an integer specifying the number of points used to construct the
//                      bounding ellipsoid.
//
// OUTPUT:
//      A double containing the value of the updated enlargement fraction.
//

double MultiEllipsoidSampler::updateEnlargementFraction(const int clusterSize)
{
    double updatedEnlargementFraction = initialEnlargementFraction * exp( shrinkingRate * logRemainingPriorMass 
                                            + 0.5 * log(static_cast<double>(NlivePoints) / clusterSize) );
    
    return updatedEnlargementFraction;
}














// MultiEllipsoidSampler::getEllipsoids()
//
// PURPOSE:
//      Gets private data member ellipsoids.
//
// OUTPUT:
//      a vector of Ellipsoids class objects containing the ellipsoids 
//      computed during the sampling process.
//

vector<Ellipsoid> MultiEllipsoidSampler::getEllipsoids()
{
    return ellipsoids;
}












// MultiEllipsoidSampler::getInitialEnlargementFraction()
//
// PURPOSE:
//      Gets private data member initialEnlargementFraction.
//
// OUTPUT:
//      a double containing the value of the initial enlargement fraction
//      of the ellipsoids.
//

double MultiEllipsoidSampler::getInitialEnlargementFraction()
{
    return initialEnlargementFraction;
}











// MultiEllipsoidSampler::getShrinkingRate()
//
// PURPOSE:
//      Gets private data member shrinkingRate.
//
// OUTPUT:
//      a double containing the value of the shrinking rate
//      of the ellipsoids.
//

double MultiEllipsoidSampler::getShrinkingRate()
{
    return shrinkingRate;
}
