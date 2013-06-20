#include "MultiEllipsoidSampler.h"

// MultiEllipsoidSampler::MultiEllipsoidSampler()
//
// PURPOSE: 
//      Class constructor.
//
// INPUT:
//      printOnTheScreen:         true if the results are to be printed on the screen, false otherwise
//      ptrPriors:                vector of Prior class objects containing the priors used in the problem.
//      likelihood:               Likelihood class object used for likelihood sampling.
//      metric:                   Metric class object to contain the metric used in the problem.
//      clusterer:                Clusterer class object specifying the type of clustering algorithm to be used.      
//      Nobjects:                 Initial number of objects coming from the main nested sampling loop
//      initialEnlargementFactor: Initial value of the enlargement for the ellipsoids
//      shrinkingRate:            Shrinking rate of the enlargement factors, between 0 and 1.
//

MultiEllipsoidSampler::MultiEllipsoidSampler(const bool printOnTheScreen, vector<Prior*> ptrPriors, 
                                             Likelihood &likelihood, Metric &metric, Clusterer &clusterer, 
                                             const int Nobjects, const double initialEnlargementFactor, const double shrinkingRate
                                            )
: NestedSampler(printOnTheScreen, Nobjects, ptrPriors, likelihood, metric, clusterer),
  initialEnlargementFactor(initialEnlargementFactor),
  shrinkingRate(shrinkingRate)
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
//      Draws an arbitrary number of points from one of the ellipsoids selected in the sample
//      for those overlapping, and starts a separate nesting iteration for each
//      ellipsoid which is isolated, i.e. non-overlapping.
//
// INPUT:
//      totalSampleOfParameters:  Eigen Array matrix of size (Ndimensions, Nobjects)
//      Nclusters:                Optimal number of clusters found by clustering algorithm
//      clusterIndices:           Indices of clusters for each point of the sample
//      logTotalWidthInPriorMass: Log Value of total prior volume from beginning to the actual nested iteration.
//      drawnSampleOfParameters:  Eigen Array matrix of size (Ndimensions,Ndraws) to contain the
//                                 coordinates of the drawn point to be used for the next nesting loop. 
//                                 When used for the first time, the array contains the coordinates of the 
//                                 worst nested object to be updated.
//      maxNdrawAttempts:         Maximum number of attempts allowed when drawing from a single ellipsoid.
//      
//
// OUTPUT:
//      void
//

void MultiEllipsoidSampler::drawWithConstraint(const RefArrayXXd totalSampleOfParameters, const int Nclusters, const vector<int> &clusterIndices, 
                                               const double logTotalWidthInPriorMass, RefArrayXXd drawnSampleOfParameters, const int maxNdrawAttempts)
{    
    assert(totalSampleOfParameters.cols() == clusterIndices.size());
    assert(drawnSampleOfParameters.rows() == totalSampleOfParameters.rows());
    assert(Nclusters > 0);

    int Ndraws = drawnSampleOfParameters.cols();


    // Compute covariance matrices and centers for ellipsoids associated with each cluster.
    // Compute eigenvalues and eigenvectors, enlarge eigenvalues and compute their hyper volumes.

    logRemainingWidthInPriorMass = log(1.0 - exp(logTotalWidthInPriorMass));
    computeEllipsoids(totalSampleOfParameters, Nclusters, clusterIndices, logRemainingWidthInPriorMass);
    
    if (printOnTheScreen)
    {   
        cerr << "=========================================" << endl;
        cerr << "Information on Ellipsoidal Sampler" << endl;
        cerr << "=========================================" << endl;
        cerr << "Nellipsoids computed: " << Nellipsoids << endl;
    }


    // Find which ellipsoids are overlapping and which are not
    
    findOverlappingEllipsoids();
    

    // Define some variables to be used in the sampling process

    bool someEllipsoidsOverlap = false;         // True if overlapping ellipsoids exist
    bool someEllipsoidsDoNotOverlap = false;    // True if non-overlapping ellipsoids exist
    bool ellipsoidIsOverlapping = false;        // True if selected ellipsoid is overlapping
    bool pointIsRejectedFromPrior;              // True if drawn point does not verify prior conditions
    bool overlapIsFound;                        // True if overlapping with another ellipsoid has been detected                   
    int ellipsoidIndex;                         // Ellipsoid index
    int ellipsoidIndex2;                        // Ellipsoid index #2
    int mergedUniformIndex;                     // Ellipsoid index among all ellipsoids (overlapping + non-overlapping)
    int NoverlappingEllipsoids;                 // Number of overlapping ellipsoids
    int NnonOverlappingEllipsoids;              // Number of non-overlapping ellipsoids
    int NdimensionsPerPrior;                    // Number of dimensions having same type of prior 
    int Nloops;                                 // Number of attempts when checking likelihood constraint
    int Noverlaps;                              // Number of detected overlaps
    int actualNdimensions;
    double volumeProbability;
    double actualProbability; 
    double rejectProbability;
    double totalVolume1 = 0;                             // Total volume of overdlapping ellipsoids
    double totalVolume2 = 0;                             // Total volume of non-overlapping ellipsoids
    double logLikelihood;
    vector<int> mergedEllipsoidsIndices;                 // Total vector of ellipsoids indices (overlapping + non-overlapping)
    ArrayXd drawnParametersPerObject(Ndimensions);       // Coordinates for the drawn point
    ArrayXd drawnParametersPerPrior;                     // Coordinates for the drawn point corresponding to one type of prior
    ArrayXd referenceParametersPerObject(Ndimensions);   // Coordinates for reference point
    ArrayXd referenceParametersPerPrior;                 // Coordinates for reference point corresponding to one type of prior
    ArrayXXd drawnAndReferenceParametersPerPrior;        // Coordinates for drawn and reference point in case of non-uniform priors   

    
    // Compute total volumes for overlapping and non-overlapping ellipsoids, if any

    if (overlappingEllipsoidsIndices[0] != -1)
    {
        NoverlappingEllipsoids = overlappingEllipsoidsIndices.size();
        
        if (printOnTheScreen)
        {    
            cerr << "Nellipsoids OVERLAPPING: " << NoverlappingEllipsoids << endl;
        }
    
        for (int m = 0; m < NoverlappingEllipsoids; m++)
        {   
            ellipsoidIndex = overlappingEllipsoidsIndices[m];
            totalVolume1 += ellipsoids[ellipsoidIndex].getHyperVolume();           // Compute total volume of overlapping ellipsoids
        }

        someEllipsoidsOverlap = true;
    }
    else
    {
        NoverlappingEllipsoids = 0;
        someEllipsoidsOverlap = false;
    }
    
    
    if (nonOverlappingEllipsoidsIndices[0] != -1)
    {
        NnonOverlappingEllipsoids = nonOverlappingEllipsoidsIndices.size();
        
        if (printOnTheScreen)
        {    
            cerr << "Nellipsoids NON-OVERLAPPING: " << NnonOverlappingEllipsoids << endl;
        }
        
        for (int m = 0; m < NnonOverlappingEllipsoids; m++)
        {   
            ellipsoidIndex = nonOverlappingEllipsoidsIndices[m];
            totalVolume2 += ellipsoids[ellipsoidIndex].getHyperVolume();           // Compute total volume of non-overlapping ellipsoids
        }

        someEllipsoidsDoNotOverlap = true;
    }
    else
    {
        NnonOverlappingEllipsoids = 0;
        someEllipsoidsDoNotOverlap = false;
    }


    // Start the do-while loop for checking maximum number of attempts when drawing from ellipsoids 

    do 
    {
        // Pick up one ellipsoid with probability according to its volume

        if (someEllipsoidsOverlap && someEllipsoidsDoNotOverlap) // If both overlapping and non-overlapping ellipsoids exist then...
        {
            // Copy the indices of the overlapping and nonoverlapping ellipsoids into 'mergedEllipsoidsIndices'

            mergedEllipsoidsIndices.resize(NoverlappingEllipsoids + NnonOverlappingEllipsoids);

            for (int i = 0; i < NoverlappingEllipsoids; ++i)
            {
                mergedEllipsoidsIndices[i] = overlappingEllipsoidsIndices[i];
            }
            
            for (int i = 0; i < NnonOverlappingEllipsoids; ++i)
            {
                mergedEllipsoidsIndices[i+NoverlappingEllipsoids] = nonOverlappingEllipsoidsIndices[i];
            }

            // Define an discrete uniform random generator

            uniform_int_distribution<int> uniformIndex(0, NoverlappingEllipsoids + NnonOverlappingEllipsoids - 1);

            do
            {
                mergedUniformIndex = uniformIndex(engine);
                ellipsoidIndex = mergedEllipsoidsIndices[mergedUniformIndex];   // Pick up one ellipsoid randomly
                
                
                // Compute probability for the selected ellipsoid
                
                volumeProbability = ellipsoids[ellipsoidIndex].getHyperVolume()/(totalVolume1+totalVolume2);
                actualProbability = uniform(engine);    // Give a probability value between 0 and 1
           

                // Determine if selected ellipsoid is among overlapping ones

                if (mergedUniformIndex < NoverlappingEllipsoids)
                {
                    ellipsoidIsOverlapping = true;
                }
                else
                {
                    ellipsoidIsOverlapping = false;
                }
            }
            while (actualProbability > volumeProbability);    // If actualProbability < volumeProbability then pick the corresponding ellipsoid
        }
        else 
            if (someEllipsoidsOverlap && !someEllipsoidsDoNotOverlap)    // If only overlapping ellipsoids exist then...
            {
                mergedEllipsoidsIndices.resize(NoverlappingEllipsoids);
                mergedEllipsoidsIndices = overlappingEllipsoidsIndices;
                uniform_int_distribution<int> uniformIndex1(0, NoverlappingEllipsoids-1);
                ellipsoidIsOverlapping = true;
        
                do
                {
                    ellipsoidIndex = mergedEllipsoidsIndices[uniformIndex1(engine)];    // Pick up one ellipsoid randomly
                    
                    
                    // Compute probability for the selected ellipsoid
                    
                    volumeProbability = ellipsoids[ellipsoidIndex].getHyperVolume()/totalVolume1;    
                    actualProbability = uniform(engine);    // Give a probability value between 0 and 1
                }
                while (actualProbability > volumeProbability);    // If actualProbability < volumeProbability then pick the corresponding ellipsoid
            }
        else 
            if (!someEllipsoidsOverlap && someEllipsoidsDoNotOverlap)    // If only non-overlapping ellipsoids exist then...
            {
                if (NnonOverlappingEllipsoids == 1)    // If only one non-overlapping ellipsoid exist, select it directly.
                    ellipsoidIndex = nonOverlappingEllipsoidsIndices[0];
                else
                {
                    mergedEllipsoidsIndices.resize(NnonOverlappingEllipsoids);
                    mergedEllipsoidsIndices = nonOverlappingEllipsoidsIndices;
                    uniform_int_distribution<int> uniformIndex2(0, NnonOverlappingEllipsoids-1);
                    ellipsoidIsOverlapping = false;
        
                    do
                    {
                        ellipsoidIndex = mergedEllipsoidsIndices[uniformIndex2(engine)];    // Pick up one ellipsoid randomly
                    
                    
                        // Compute probability for the selected ellipsoid
                    
                        volumeProbability = ellipsoids[ellipsoidIndex].getHyperVolume()/totalVolume2;    
                        actualProbability = uniform(engine);    // Give a probability value between 0 and 1
                    }
                    while (actualProbability > volumeProbability);    // If actualProbability < volumeProbability then pick the corresponding ellipsoid
                }
            }

        actualNdimensions = 0;
        Nloops = 0;

        for (int m = 0; m < Ndraws; m++)    // FOR loop over the total number of points to be drawn from one ellipsoid - Default is Ndraws = 1
        {
            drawnParametersPerObject = drawnSampleOfParameters.col(m);
            
            if (!ellipsoidIsOverlapping)    // Isolated ellipsoid sampling
            {
                // Start do-while loop for checking likelihood constraint when drawing from ellipsoids

                do   
                {
                    // Start do-while loop for checking prior distribution when drawing from ellipsoids

                    do 
                    {
                        // Draw one point from the chosen ellipsoid

                        ellipsoids[ellipsoidIndex].drawPoint(drawnParametersPerObject);
                        pointIsRejectedFromPrior = false;
                      

                        // Split the sample of drawn parameters according to type of priors and check if prior
                        // conditions are verified

                        for (int i = 0; i < ptrPriors.size(); i++)
                        {
                            NdimensionsPerPrior = ptrPriors[i]->getNdimensions();
                            drawnParametersPerPrior.resize(NdimensionsPerPrior);
                            drawnParametersPerPrior = drawnParametersPerObject.segment(actualNdimensions,NdimensionsPerPrior);      
                            

                            // Check if the prior type selected is uniform. In case it is not, draw a second
                            // point from the ellipsoid for accomplishing the sampling from the prior.

                            if (ptrPriors[i]->priorIsUniform())
                            {
                                // Only for uniform priors

                                pointIsRejectedFromPrior += ptrPriors[i]->pointIsRejected(drawnParametersPerPrior);
                            }
                            else 
                                if (!(ptrPriors[i]->priorIsUniform()))
                                {
                                    // Only for non-uniform priors

                                    ellipsoids[ellipsoidIndex].drawPoint(referenceParametersPerObject);
                                    drawnAndReferenceParametersPerPrior.resize(NdimensionsPerPrior,2);
                                    referenceParametersPerPrior = referenceParametersPerObject.segment(actualNdimensions,NdimensionsPerPrior);      
                                    drawnAndReferenceParametersPerPrior.col(0) = drawnParametersPerPrior;      
                                    drawnAndReferenceParametersPerPrior.col(1) = referenceParametersPerPrior;
                                    pointIsRejectedFromPrior += ptrPriors[i]->pointIsRejected(drawnAndReferenceParametersPerPrior);
                                }
                            
                            actualNdimensions += NdimensionsPerPrior;       // Move to next prior dimensions
                        }


                        // If sampling is not verified for at least one type of prior, repeat the drawing of the new object 
                        
                        actualNdimensions = 0;
                        pointIsRejectedFromPrior = static_cast<bool>(pointIsRejectedFromPrior);
                    } 
                    while (pointIsRejectedFromPrior);        // Rejection according to prior distribution
                
                    //cout << "From isolated: " << Nloops << endl;
                    logLikelihood = likelihood.logValue(drawnParametersPerObject);
                    Nloops++;     
                }
                while ((logLikelihood <= actualLogLikelihoodConstraint) && (Nloops <= maxNdrawAttempts)); // Rejection according to likelihood constraint
        
            }   // End sampling from isolated ellipsoid
            else 
                if (ellipsoidIsOverlapping)        // Overlapping ellipsoid sampling
                {
                    // Sample uniformly from chosen overlapping ellipsoid until new parameter with 
                    // logLikelihood > actualLogLikelihoodConstraint is found.
                    // If parameter is contained in Noverlaps ellipsoids, then accept parameter with probability 1/Noverlaps.
        
                    // Start do-while loop for sampling from overlapping ellipsoids 

                    do  
                    {
                        // Start do-while loop for checking likelihood constraint when drawing from ellipsoids

                        do  
                        {
                            // Start do-while loop for checking prior distribution when drawing from ellipsoids

                            do 
                            {
                                // Draw one point from the chosen ellipsoid

                                ellipsoids[ellipsoidIndex].drawPoint(drawnParametersPerObject);
                                pointIsRejectedFromPrior = false;
                      

                                // Split the sample of drawn parameters according to type of priors and check if prior
                                // conditions are verified

                                for (int i = 0; i < ptrPriors.size(); i++)
                                {
                                    NdimensionsPerPrior = ptrPriors[i]->getNdimensions();
                                    drawnParametersPerPrior.resize(NdimensionsPerPrior);
                                    drawnParametersPerPrior = drawnParametersPerObject.segment(actualNdimensions,NdimensionsPerPrior);      
                            

                                    // Check if the prior type selected is uniform. In case it is not, draw a second
                                    // point from the ellipsoid for accomplishing the sampling from the prior.

                                    if (ptrPriors[i]->priorIsUniform())
                                    {
                                        // Only for uniform priors

                                        pointIsRejectedFromPrior += ptrPriors[i]->pointIsRejected(drawnParametersPerPrior);
                                    }
                                    else 
                                        if (!(ptrPriors[i]->priorIsUniform()))
                                        {
                                            // Only for non-uniform priors

                                            ellipsoids[ellipsoidIndex].drawPoint(referenceParametersPerObject);
                                            drawnAndReferenceParametersPerPrior.resize(NdimensionsPerPrior,2);
                                            referenceParametersPerPrior = referenceParametersPerObject.segment(actualNdimensions,NdimensionsPerPrior);
                                            drawnAndReferenceParametersPerPrior.col(0) = drawnParametersPerPrior;      
                                            drawnAndReferenceParametersPerPrior.col(1) = referenceParametersPerPrior;
                                            pointIsRejectedFromPrior += ptrPriors[i]->pointIsRejected(drawnAndReferenceParametersPerPrior);
                                        }
                            
                                    actualNdimensions += NdimensionsPerPrior;       // Move to next prior dimensions
                                }
                    

                                // If sampling is not verified for at least one type of prior, repeat the drawing of the new object 
                        
                                actualNdimensions = 0;
                                pointIsRejectedFromPrior = static_cast<bool>(pointIsRejectedFromPrior);
                            } 
                            while (pointIsRejectedFromPrior);        // Rejection according to prior distribution
                         
                            //cout << "From overlapping: " << Nloops << endl;
                            logLikelihood = likelihood.logValue(drawnParametersPerObject);
                            Nloops++;    
                    
                        }
                        while ((logLikelihood <= actualLogLikelihoodConstraint) && (Nloops <= maxNdrawAttempts));         
                
                        if ((Nloops >= maxNdrawAttempts) && (logLikelihood <= actualLogLikelihoodConstraint))
                            break;
                
                        Noverlaps = 1;        // Start with self-overlap only
                        overlapIsFound = false;

                        for (int i = 0; i < NoverlappingEllipsoids; i++)        // Check other possibile overlapping ellipsoids
                        {
                            if (overlappingEllipsoidsIndices[i] == ellipsoidIndex)     // Skip if self-overlap
                                continue;
                            else
                            {
                                // Check if point belongs to ellipsoid #2

                                ellipsoidIndex2 = overlappingEllipsoidsIndices[i];

                                if (ellipsoids[ellipsoidIndex2].containsPoint(drawnParametersPerObject))
                                    Noverlaps++;
                                else
                                    continue;
                            }
                        }

                        if (Noverlaps == 1)                             // If no overlaps found, go to next point m
                            break;
                        else
                        {
                            rejectProbability = 1./Noverlaps;           // Set rejection probability value
                            actualProbability = uniform(engine);        // Evaluate actual probability value
                        }
                
                    }
                    while (actualProbability > rejectProbability);      // If actual probability value < rejection probability then...
                                                                        // ...accept point and end function
                } // End sampling from overlapping ellipsoids
 
            if ((Nloops >= maxNdrawAttempts) && (logLikelihood <= actualLogLikelihoodConstraint))
            {
                if ((NnonOverlappingEllipsoids == 1) && (NoverlappingEllipsoids == 0))
                {
                    cerr << "Drawing from single ellipsoid stuck into flat likelihood region." << endl; 
                    cerr << "Quitting program." << endl;
                    exit(EXIT_FAILURE);
                }
                else
                    break;
            }
       
            drawnSampleOfParameters.col(m) = drawnParametersPerObject;  // Update set of parameters for one object into total sample
        } // End FOR loop
    }
    while ((Nloops >= maxNdrawAttempts) && (logLikelihood <= actualLogLikelihoodConstraint));     // Restart selection of ellipsoid in case...
                                                                                                  // ...the number of attempts (Nloops) exceeds the limit
}












// MultiEllipsoidSampler::computeEllipsoids()
//
// PURPOSE:
//      Computes covariance matrices and center coordinates of all the ellipsoids
//      associated to each cluster of the sample. The egeivalues decomposition
//      is also done and eigenvalues are enlarged afterwards. All the results
//      are stored in the private data members.
//
// INPUT:
//      totalSampleOfParameters:      Eigen Array matrix of size (Ndimensions, Nobjects)
//                                     containing the total sample of points to be split into clusters.
//      clusterIndices:               One dimensional Eigen Array containing the integer indices
//                                     of the clusters as obtained from the clustering algorithm.
//      logRemainingWidthInPriorMass: Log Value of the remaining prior volume at the actual nested iteration.
//
// OUTPUT:
//      void
//

void MultiEllipsoidSampler::computeEllipsoids(const RefArrayXXd totalSampleOfParameters, const int Nclusters, const vector<int> &clusterIndices, const double logRemainingWidthInPriorMass)
{
    assert(totalSampleOfParameters.cols() == clusterIndices.size());
    assert(totalSampleOfParameters.cols() > Ndimensions + 1);            // At least Ndimensions + 1 points are required.
 
    NobjectsPerCluster.assign(Nclusters, 0);

    // Divide the sample according to the clustering done

    // Find number of points (objects) per cluster

    for (int i = 0; i < Nobjects; i++) 
    {
        NobjectsPerCluster[clusterIndices[i]]++;
    }

    ArrayXd oneDimensionSampleOfParameters(Nobjects);
    ArrayXXd totalSampleOfParametersOrdered = totalSampleOfParameters;
    vector<int> clusterIndicesCopy(clusterIndices);

    // Order points in each dimension according to increasing cluster indices

    for (int i = 0; i < Ndimensions; i++)     
    {
        oneDimensionSampleOfParameters = totalSampleOfParameters.row(i);
        Functions::sortElementsInt(clusterIndicesCopy, oneDimensionSampleOfParameters);
        totalSampleOfParametersOrdered.row(i) = oneDimensionSampleOfParameters;
        clusterIndicesCopy = clusterIndices;
    }


    // Build vector of ellipsoids objects for all the clusters found in the total sample having 
    // a number of points not lower than Ndimensions + 1.
    // In each object, covariance matrix, center coordinates, eigenvalues and eigenvectoes are computed,
    // according to the corresponding enlargement factor.
    
    ArrayXXd clusterSample;
    int actualNobjects = 0;
    double logEnlargementFactor;
    double enlargementFactor;
    Nellipsoids = 0;                    // The total number of ellipsoids computed

    for (int i = 0; i < Nclusters; i++)
    {   
        // Skip cluster if number of points is not large enough

        if (NobjectsPerCluster[i] <= Ndimensions + 1) 
        {
            actualNobjects += NobjectsPerCluster[i];
            continue;
        }
        else
        {
            clusterSample.resize(Ndimensions, NobjectsPerCluster[i]);
            clusterSample = totalSampleOfParametersOrdered.block(0, actualNobjects, Ndimensions, NobjectsPerCluster[i]);
            actualNobjects += NobjectsPerCluster[i];

            // Compute the enlargement factor

            logEnlargementFactor = log(initialEnlargementFactor) + shrinkingRate * logRemainingWidthInPriorMass 
                                                                 + 0.5 * log(Nobjects/NobjectsPerCluster[i]);
            enlargementFactor = exp(logEnlargementFactor);

            // Insert ellipsoid in our vector

            ellipsoids.insert(ellipsoids.end(), Ellipsoid(clusterSample, enlargementFactor));
            Nellipsoids++;
        }
    }
}













// MultiEllipsoidSampler::findOverlappingEllipsoids()
//
// PURPOSE:
//      Determines which of Nclusters input ellipsoids are overlapping and which
//      are not. The results are stored in the private data members specifying the 
//      subscripts of the ellipsoids that are not overlapping. 
//      A single element array containing the value -1 is saved 
//      in case no matchings are found.
//
// OUTPUT:
//      An Eigen Array of integers specifying the subscripts of the ellipsoids that 
//      are not overlapping is initialized as protected member. 
//      A single element array containing the value -1 is returned in case no 
//      non-overlapping ellipsoids are found.
//
void MultiEllipsoidSampler::findOverlappingEllipsoids()
{
    overlappingEllipsoidsIndices.assign(1, -1);       // -1: Start with no overlapping ellipsoids found
    nonOverlappingEllipsoidsIndices.assign(1, -1);    // -1: Start with no non-overlapping ellipsoids found

    // If only one ellipsoid is found, then return one non-overlapping ellipsoid

    if (Nellipsoids == 1) 
    {
        nonOverlappingEllipsoidsIndices[0] = 0;
        return;
    }

    bool saveFlagI = true;
    bool saveFlagJ = true;
    int countOverlap;
    int countIndex = 0;
    

    // Identify the overlapping ellipsoids

    for (int i = 0; i < Nellipsoids-1; i++)
    {
        saveFlagI = true;       // Reset save index flag for i
        countOverlap = 0;       // Reset overlap counter
        
        for (int j = i + 1; j < Nellipsoids; j++)
        {   
            saveFlagJ = true;        // Reset save index flag for j

            if (ellipsoids[i].overlapsWith(ellipsoids[j]))  
            {
                countOverlap++;

                for (int k = 0; k < overlappingEllipsoidsIndices.size(); k++) // Check if j is saved already
                {
                    if (j == overlappingEllipsoidsIndices[k])       
                    {
                        saveFlagJ = false; // If j is saved already, don't save it again and go to next j
                        break;
                    }
                }

                
                if (saveFlagJ)      // If j is not saved yet, then save it
                {
                    overlappingEllipsoidsIndices.resize(countIndex+1);
                    overlappingEllipsoidsIndices[countIndex] = j;
                    countIndex++;
                }
            }
        }

        if (i == 0 && countOverlap != 0)   // If first i and at least one overlap is found, save also i and go to next i
        {
            overlappingEllipsoidsIndices.resize(countIndex+1);
            overlappingEllipsoidsIndices[countIndex] = i;
            countIndex++;
            continue;
        }
        else 
            if (i > 0 && countOverlap != 0)     // If i is not the first one and at least one overlap is found ...
            {
                for (int k = 0; k < overlappingEllipsoidsIndices.size(); k++) // Check if i is saved already
                {
                    if (i == overlappingEllipsoidsIndices[k])       
                    {
                        saveFlagI = false;       // If i is saved already, don't save it again and go to next i
                        break;
                    }
                }

                if (saveFlagI)      // If i is not saved yet, then save it
                {
                    overlappingEllipsoidsIndices.resize(countIndex+1);
                    overlappingEllipsoidsIndices[countIndex] = i;
                    countIndex++;
                }
            }
    }

    int Noverlaps;

    if (overlappingEllipsoidsIndices[0] != -1)
        Noverlaps = overlappingEllipsoidsIndices.size();
    else
        Noverlaps = 0;

    int NnonOverlaps = Nellipsoids - Noverlaps;

    if (NnonOverlaps == 0)          // If no non-overlapping ellipsoids are found, keep vector element to -1
        return;
    else
    {
        nonOverlappingEllipsoidsIndices.resize(NnonOverlaps);       // Otherwise, at least 1 NnonOverlaps ellipsoid is found
        countIndex = 0;

        for (int i = 0; i < Nellipsoids; i++)
        {
            saveFlagI = true;       // Reset save flag for i

            for (int j = 0; j < Noverlaps; j++)     // Check which ellipsoids are already overlapping
            {
                if (i == overlappingEllipsoidsIndices[j])       // If ellipsoid i-th already overlap, then go to next i-th ellipsoid
                {
                    saveFlagI = false;
                    break;
                }
            }

            if (saveFlagI)      // If ellipsoid is not overlapping, save it among non-overlapping ellipsoids
            {   
                nonOverlappingEllipsoidsIndices[countIndex] = i;
                countIndex++;
            }
        }
    }

}









// MultiEllipsoidSampler::getNonOverlappingEllipsoidsIndices()
//
// PURPOSE:
//      Gets private data member nonOverlappingEllipsoidsIndices.
//
// OUTPUT:
//      an Eigen Array containing the indices of the ellipsoids that
//      are not overlapping.
//

vector<int> MultiEllipsoidSampler::getNonOverlappingEllipsoidsIndices()
{
    return nonOverlappingEllipsoidsIndices;
}














// MultiEllipsoidSampler::getOverlappingEllipsoidsIndices()
//
// PURPOSE:
//      Gets private data member overlappingEllipsoidsIndices.
//
// OUTPUT:
//      an Eigen Array containing the indices of the ellipsoids that
//      are overlapping.
//

vector<int> MultiEllipsoidSampler::getOverlappingEllipsoidsIndices()
{
    return overlappingEllipsoidsIndices;
}














// MultiEllipsoidSampler::getEllipsoidsVector()
//
// PURPOSE:
//      Gets private data member ellipsoids.
//
// OUTPUT:
//      a vector of Ellipsoids class objects containing the ellipsoids 
//      computed during the sampling process.
//

vector<Ellipsoid> MultiEllipsoidSampler::getEllipsoidsVector()
{
    return ellipsoids;
}
