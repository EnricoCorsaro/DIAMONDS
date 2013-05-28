#include "MultiEllipsoidSampler.h"

// MultiEllipsoidSampler::MultiEllipsoidSampler()
//
// PURPOSE: 
//      Base class constructor.
//
// INPUT:
//      ptrPriorsVector: a vector of Prior class objects containing the priors used in the problem.
//      likelihood: a Likelihood class object used for likelihood sampling.
//      metric: a Metric class object to contain the metric used in the problem.
//      clusterer: a Clusterer class object specifying the type of clustering algorithm to be used.      
//      Nobjects: the initial number of objects coming from the main nested sampling loop
//      initialEnlargementFactor: initial value of the enlargement for the ellipsoids
//      alpha: a value between 0 and 1 that defines the shrinking rate of the enlargement factors
//

MultiEllipsoidSampler::MultiEllipsoidSampler(vector<Prior*> ptrPriorsVector, Likelihood &likelihood, Metric &metric, Clusterer &clusterer, const int Nobjects, const double initialEnlargementFactor, const double alpha)
: NestedSampler(Nobjects, ptrPriorsVector, likelihood, metric, clusterer),
  uniform(0.0,1.0),
  normal(0.0,1.0),
  initialEnlargementFactor(initialEnlargementFactor),
  alpha(alpha)
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
//      Draws an arbitrary number of objects from one of the ellipsoids selected in the sample
//      for those overlapping, and starts a separate nesting iteration for each
//      ellipsoid which is isolated, i.e. non-overlapping.
//
// INPUT:
//      totalSampleOfParameters: an Eigen Array matrix of size (Ndimensions, Nobjects)
//      Nclusters: optimal number of clusters found by clustering algorithm
//      clusterIndices: indices of clusters for each point of the sample
//      logTotalWidthInPriorMass: log Value of total prior volume from beginning to the actual nested iteration.
//      drawnSampleOfParameters: an Eigen Array matrix of size (Ndimensions,Ndraws) to contain the
//      coordinates of the drawn point to be used for the next nesting loop. When used
//      for the first time, the array contains the coordinates of the worst nested
//      object to be updated.
//      
//
// OUTPUT:
//      void
//

void MultiEllipsoidSampler::drawWithConstraint(const RefArrayXXd totalSampleOfParameters, const int Nclusters, const RefArrayXi clusterIndices, const double logTotalWidthInPriorMass, RefArrayXXd drawnSampleOfParameters)
{    
    assert(totalSampleOfParameters.cols() == clusterIndices.size());
    assert(drawnSampleOfParameters.rows() == totalSampleOfParameters.rows());
    assert(Nclusters > 0);

    int Ndraws = drawnSampleOfParameters.cols();


    // Compute covariance matrices and centers for ellipsoids associated with each cluster.
    // Compute eigenvalues and eigenvectors, enlarge eigenvalues and compute their hyper volumes.

    logRemainingWidthInPriorMass = log(1.0 - exp(logTotalWidthInPriorMass));
    computeEllipsoids(totalSampleOfParameters, Nclusters, clusterIndices, logRemainingWidthInPriorMass);
    
    
    // Find which ellipsoids are overlapping and which are not
    
    findOverlappingEllipsoids();
    

    // Define some variables to be used in the sampling process

    bool ellipsoidsOverlap; 
    bool ellipsoidsNonOverlap;
    bool ellipsoidIsInOverlap;
    bool pointIsRejectedFromPrior;
    int ellipsoidIndex;                   // Ellipsoid index
    int ellipsoidIndex2;
    int NoverlappingEllipsoids;
    int NnonOverlappingEllipsoids;
    int NdimensionsPerPrior;                // Number of dimensions having same type of prior 
    int Nloops;
    int NloopMaximum = 200;              // Maximum number of iterations for likelihood constraint
    int actualNdimensions;
    double volumeProbability;
    double actualProbability; 
    double rejectProbability;
    double totalVolume1;
    double totalVolume2;
    double logLikelihood;
    ArrayXi mergedEllipsoidsIndices;
    ArrayXd drawnParametersPerObject(Ndimensions);       // Coordinates for the drawn point
    ArrayXd drawnParametersPerPrior;                     // Coordinates for the drawn point corresponding to one type of prior
    ArrayXd referenceParametersPerObject(Ndimensions);   // Coordinates for reference point
    ArrayXd referenceParametersPerPrior;                 // Coordinates for reference point corresponding to one type of prior
    ArrayXXd drawnAndReferenceParametersPerPrior;        // Coordinates for drawn and reference point in case of non-uniform priors   

    ellipsoidsOverlap = false;
    ellipsoidsNonOverlap = false;
    ellipsoidIsInOverlap = false;
    totalVolume1 = 0;        // For overlapping ellipsoids
    totalVolume2 = 0;        // For non-overlapping ellipsoids
    
    
    // Choose one ellipsoid with probability according to its volume fraction

    if (overlappingEllipsoidsIndices(0) != -1)
    {
        NoverlappingEllipsoids = overlappingEllipsoidsIndices.size();
    
        for (int m = 0; m < NoverlappingEllipsoids; m++)
        {   
            ellipsoidIndex = overlappingEllipsoidsIndices(m);
            totalVolume1 += ellipsoidsVector[ellipsoidIndex].getHyperVolume();           // Compute total volume of overlapping ellipsoids
        }

        ellipsoidsOverlap = true;
    }
    else
    {
        NoverlappingEllipsoids = 0;
        ellipsoidsOverlap = false;
    }
    
    
    if (nonOverlappingEllipsoidsIndices(0) != -1)
    {
        NnonOverlappingEllipsoids = nonOverlappingEllipsoidsIndices.size();
        
        for (int m = 0; m < NnonOverlappingEllipsoids; m++)
        {   
            ellipsoidIndex = nonOverlappingEllipsoidsIndices(m);
            totalVolume2 += ellipsoidsVector[ellipsoidIndex].getHyperVolume();           // Compute total volume of non-overlapping ellipsoids
        }

        ellipsoidsNonOverlap = true;
    }
    else
    {
        NnonOverlappingEllipsoids = 0;
        ellipsoidsNonOverlap = false;
    }

   
    do 
    {
        // Pick up one ellipsoid with probability according to its volume

        if (ellipsoidsOverlap && ellipsoidsNonOverlap)
        {
            mergedEllipsoidsIndices.resize(NoverlappingEllipsoids+NnonOverlappingEllipsoids);
            mergedEllipsoidsIndices.head(NoverlappingEllipsoids) = overlappingEllipsoidsIndices;
            mergedEllipsoidsIndices.segment(NoverlappingEllipsoids, NnonOverlappingEllipsoids) = nonOverlappingEllipsoidsIndices;
            uniform_int_distribution<int> uniformIndex(0, NoverlappingEllipsoids+NnonOverlappingEllipsoids-1);

            int mergedUniformIndex;

            do
            {
                mergedUniformIndex = uniformIndex(engine);
                ellipsoidIndex = mergedEllipsoidsIndices(mergedUniformIndex);      // Pick up one ellipsoid randomly
                volumeProbability = ellipsoidsVector[ellipsoidIndex].getHyperVolume()/(totalVolume1+totalVolume2);// Compute probability for the selected ellipsoid
                actualProbability = uniform(engine);    // Give a probability value between 0 and 1
            
                if (mergedUniformIndex < NoverlappingEllipsoids)
                {
                    ellipsoidIsInOverlap = true;
                }
                else
                {
                    ellipsoidIsInOverlap = false;
                }
            }
            while (actualProbability > volumeProbability);      // If actualProbability < volumeProbability then pick the corresponding ellipsoid
        }
        else if (ellipsoidsOverlap && !ellipsoidsNonOverlap)
        {
            mergedEllipsoidsIndices.resize(NoverlappingEllipsoids);
            mergedEllipsoidsIndices = overlappingEllipsoidsIndices;
            uniform_int_distribution<int> uniformIndex1(0, NoverlappingEllipsoids-1);
            ellipsoidIsInOverlap = true;
        
            do
            {
                ellipsoidIndex = mergedEllipsoidsIndices(uniformIndex1(engine));      // Pick up one ellipsoid randomly
                volumeProbability = ellipsoidsVector[ellipsoidIndex].getHyperVolume()/totalVolume1;    // Compute probability for the selected ellipsoid
                actualProbability = uniform(engine);                                // Give a probability value between 0 and 1
            }
            while (actualProbability > volumeProbability);      // If actualProbability < volumeProbability then pick the corresponding ellipsoid
        }
        else if (!ellipsoidsOverlap && ellipsoidsNonOverlap)
        {
            mergedEllipsoidsIndices.resize(NnonOverlappingEllipsoids);
            mergedEllipsoidsIndices = nonOverlappingEllipsoidsIndices;
            uniform_int_distribution<int> uniformIndex2(0, NnonOverlappingEllipsoids-1);
            ellipsoidIsInOverlap = false;
        
            do
            {
                ellipsoidIndex = mergedEllipsoidsIndices(uniformIndex2(engine));      // Pick up one ellipsoid randomly
                volumeProbability = ellipsoidsVector[ellipsoidIndex].getHyperVolume()/totalVolume2;    // Compute probability for the selected ellipsoid
                actualProbability = uniform(engine);                                // Give a probability value between 0 and 1
            }
            while (actualProbability > volumeProbability);      // If actualProbability < volumeProbability then pick the corresponding ellipsoid
        }

        actualNdimensions = 0;
        Nloops = 0;

        if (!ellipsoidIsInOverlap)    // Isolated ellipsoid sampling
        {
            /* START OF TO BE FIXED IN FUTURE 
            // If non-overlapping ellipsoids exist start separate nested sampling 
            // from each non overlapping enlarged ellipsoid 
    
            int NobjectsInEllipsoid;              // Number of objects for the nested sampling within the identified ellipsoid
            double logRemainingWidthInPriorMassReduced;      // Shrinked prior volume for separate nesting of ellipsoid
            double logTotalEvidenceFromIsolatedEllipsoids = -DBL_MAX;  // The log value of the total contribute of local evidences from isolated ellipsoids
            double logLocalEvidence;
            END OF TO BE FIXED IN FUTURE */ 
        
            for (int m = 0; m < Ndraws; m++)
            {
                drawnParametersPerObject = drawnSampleOfParameters.col(m);
            
                do
                {
                    do
                    {
                        /* START OF TO BE FIXED IN FUTURE 
                        NobjectsInEllipsoid = ellipsoidsVector[ellipsoidIndex].getNobjects();       // Find the number of live points inside the chosen ellipsoid
                        logRemainingWidthInPriorMassReduced = logRemainingWidthInPriorMass + log(NobjectsInEllipsoid/Nobjects); // Set reduced remaining prior mass
                        cout << "Start separate nesting loop for ellipsoid #" << ellipsoidIndex << endl;


                        // Perform separate nesting iteration for the actual ellipsoid

                        logLocalEvidence = -DBL_MAX;            // Set initial local evidence to zero
                        nestFromIsolatedEllipsoid(ellipsoidIndex, logLocalEvidence);
                        logTotalEvidenceFromIsolatedEllipsoids = Functions::logExpSum(logTotalEvidenceFromIsolatedEllipsoids,logLocalEvidence)
                        END OF TO BE FIXED IN FUTURE */
                

                        // Draw one point from the chosen ellipsoid

                        drawFromHyperSphere(ellipsoidsVector[ellipsoidIndex], drawnParametersPerObject);
                        pointIsRejectedFromPrior = false;
                      

                        // Split the sample of drawn parameters according to type of priors and check if prior
                        // conditions are verified

                        for (int i = 0; i < ptrPriorsVector.size(); i++)
                        {
                            NdimensionsPerPrior = ptrPriorsVector[i]->getNdimensions();
                            drawnParametersPerPrior.resize(NdimensionsPerPrior);
                            drawnParametersPerPrior = drawnParametersPerObject.segment(actualNdimensions,NdimensionsPerPrior);      
                            

                            // Check if the prior type selected is uniform. In case it is not, draw a second
                            // point from the ellipsoid for accomplishing the sampling from the prior.

                            if (ptrPriorsVector[i]->priorIsUniform())
                            {
                                // Only for uniform priors

                                pointIsRejectedFromPrior += ptrPriorsVector[i]->pointIsRejected(drawnParametersPerPrior);
                            }
                            else
                            {
                                // Only for non-uniform priors

                                drawFromHyperSphere(ellipsoidsVector[ellipsoidIndex], referenceParametersPerObject);
                                drawnAndReferenceParametersPerPrior.resize(NdimensionsPerPrior,2);
                                referenceParametersPerPrior = referenceParametersPerObject.segment(actualNdimensions,NdimensionsPerPrior);      
                                drawnAndReferenceParametersPerPrior.col(0) = drawnParametersPerPrior;      
                                drawnAndReferenceParametersPerPrior.col(1) = referenceParametersPerPrior;
                                pointIsRejectedFromPrior += ptrPriorsVector[i]->pointIsRejected(drawnAndReferenceParametersPerPrior);
                            }
                            
                            actualNdimensions += NdimensionsPerPrior;       // Move to next prior dimensions
                        }

                        // If sampling is not verified for at least one type of prior, repeat the drawing of the new object 
                        
                        actualNdimensions = 0;
                        pointIsRejectedFromPrior = static_cast<bool>(pointIsRejectedFromPrior);
                    } 
                    while (pointIsRejectedFromPrior);        // Rejection according to prior distribution
                
                    Nloops++;
        
                    if (Nloops >= NloopMaximum)
                        break;
                
                    logLikelihood = likelihood.logValue(drawnParametersPerObject);
                }
                while (logLikelihood <= actualLogLikelihoodConstraint);         // Rejection according to likelihood constraint

                if (Nloops >= NloopMaximum)
                    break;
            
                drawnSampleOfParameters.col(m) = drawnParametersPerObject;  // Update set of parameters for one object into total sample
            }
        }
        else if (ellipsoidIsInOverlap)        // Overlapping ellipsoid sampling
        {
            // Create arrays for a second ellipsoid to check for overlap

            actualNdimensions = 0;


            // Sample uniformly from chosen ellipsoid until new parameter with logLikelihood > actualLogLikelihoodConstraint is found
            // If parameter is contained in Noverlaps ellipsoids, then accept parameter with probability 1/Noverlaps
        
            for (int m = 0; m < Ndraws; m++)
            {
                drawnParametersPerObject = drawnSampleOfParameters.col(m);
            
                do
                {
                    do
                    {
                        do
                        {
                            // Draw one point from the chosen ellipsoid

                            drawFromHyperSphere(ellipsoidsVector[ellipsoidIndex], drawnParametersPerObject);
                            pointIsRejectedFromPrior = false;
                      

                            // Split the sample of drawn parameters according to type of priors and check if prior
                            // conditions are verified

                            for (int i = 0; i < ptrPriorsVector.size(); i++)
                            {
                                NdimensionsPerPrior = ptrPriorsVector[i]->getNdimensions();
                                drawnParametersPerPrior.resize(NdimensionsPerPrior);
                                drawnParametersPerPrior = drawnParametersPerObject.segment(actualNdimensions,NdimensionsPerPrior);      
                            

                                // Check if the prior type selected is uniform. In case it is not, draw a second
                                // point from the ellipsoid for accomplishing the sampling from the prior.

                                if (ptrPriorsVector[i]->priorIsUniform())
                                {
                                    // Only for uniform priors

                                    pointIsRejectedFromPrior += ptrPriorsVector[i]->pointIsRejected(drawnParametersPerPrior);
                                }
                                else
                                {
                                    // Only for non-uniform priors

                                    drawFromHyperSphere(ellipsoidsVector[ellipsoidIndex], referenceParametersPerObject);
                                    drawnAndReferenceParametersPerPrior.resize(NdimensionsPerPrior,2);
                                    referenceParametersPerPrior = referenceParametersPerObject.segment(actualNdimensions,NdimensionsPerPrior);      
                                    drawnAndReferenceParametersPerPrior.col(0) = drawnParametersPerPrior;      
                                    drawnAndReferenceParametersPerPrior.col(1) = referenceParametersPerPrior;
                                    pointIsRejectedFromPrior += ptrPriorsVector[i]->pointIsRejected(drawnAndReferenceParametersPerPrior);
                                }
                            
                                actualNdimensions += NdimensionsPerPrior;       // Move to next prior dimensions
                            }
                    

                            // If sampling is not verified for at least one type of prior, repeat the drawing of the new object 
                        
                            actualNdimensions = 0;
                            pointIsRejectedFromPrior = static_cast<bool>(pointIsRejectedFromPrior);
                        } 
                        while (pointIsRejectedFromPrior);        // Rejection according to prior distribution
                   
                        Nloops++;
                    
                        if (Nloops >= NloopMaximum)
                            break;
                    
                        logLikelihood = likelihood.logValue(drawnParametersPerObject);
                    }
                    while (logLikelihood <= actualLogLikelihoodConstraint);         // Rejection according to likelihood constraint
                
                    if (Nloops >= NloopMaximum)
                        break;
                
                    int Noverlaps = 1;        // Start with self-overlap only
                    bool overlapFlag = false;

                    for (int i = 0; i < NoverlappingEllipsoids; i++)        // Check other possibile overlapping ellipsoids
                    {
                        if (overlappingEllipsoidsIndices(i) == ellipsoidIndex)     // Skip if self-overlap
                            continue;
                        else
                        {
                            // Check if point belongs to ellipsoid 2

                            ellipsoidIndex2 = overlappingEllipsoidsIndices(i);
                            overlapFlag = pointIsInOverlap(ellipsoidsVector[ellipsoidIndex2], drawnParametersPerObject);
                
                            if (overlapFlag)
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
                while (actualProbability > rejectProbability);      // If actual probability value < rejection probability then accept point and end function
 
                if (Nloops >= NloopMaximum)
                    break;
       
                drawnSampleOfParameters.col(m) = drawnParametersPerObject;  // Update set of parameters for one object into total sample
            }
        } // END sampling from overlapping ellipsoids
    }
    while (Nloops >= NloopMaximum);     // Select ellipsoid again if Nloops exceeds the limit

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
//      totalSampleOfParameters: an Eigen Array matrix of size (Ndimensions, Nobjects)
//      containing the total sample of points to be split into clusters.
//      clusterIndices: a one dimensional Eigen Array containing the integers
//      indices of the clusters as obtained from the clustering algorithm.
//      logRemainingWidthInPriorMass: log Value of the remaining prior volume at the actual nested iteration.
//
// OUTPUT:
//      void
//

void MultiEllipsoidSampler::computeEllipsoids(const RefArrayXXd totalSampleOfParameters, const int Nclusters, const RefArrayXi clusterIndices, const double logRemainingWidthInPriorMass)
{
    assert(totalSampleOfParameters.cols() == clusterIndices.size());
    assert(totalSampleOfParameters.cols() > Ndimensions + 1);            // At least Ndimensions + 1 points are required.
    NobjectsPerCluster.resize(Nclusters);
    NobjectsPerCluster = ArrayXi::Zero(Nclusters);                      // Fundamental initialization to zero!

    // Divide the sample according to the clustering done

    for (int i = 0; i < Nobjects; i++)         // Find number of points (objects) per cluster
    {
        NobjectsPerCluster(clusterIndices(i))++;
    }

    ArrayXd oneDimensionSampleOfParameters(Nobjects);
    ArrayXXd totalSampleOfParametersOrdered = totalSampleOfParameters;
    ArrayXi clusterIndicesCopy(Nobjects);
    clusterIndicesCopy = clusterIndices;

    for (int i = 0; i < Ndimensions; i++)     // Order points in each dimension according to increasing cluster indices
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
        if (NobjectsPerCluster(i) <= Ndimensions + 1)        // Skip cluster if number of points is not large enough
        {
            actualNobjects += NobjectsPerCluster(i);
            continue;
        }
        else
        {
            clusterSample.resize(Ndimensions, NobjectsPerCluster(i));
            clusterSample = totalSampleOfParametersOrdered.block(0, actualNobjects, Ndimensions, NobjectsPerCluster(i));
            actualNobjects += NobjectsPerCluster(i);

            Ellipsoid ellipsoid(clusterSample, i);              // Initialize ellipsoid object and save original cluster index
            ellipsoidsVector.insert(ellipsoidsVector.end(),ellipsoid);           // Resize ellipsoids vector to the actual number of ellipsoids
            logEnlargementFactor = log(initialEnlargementFactor) + alpha*logRemainingWidthInPriorMass + 0.5*log(Nobjects/NobjectsPerCluster(i));
            enlargementFactor = exp(logEnlargementFactor);
            ellipsoidsVector[Nellipsoids].build(enlargementFactor);
            Nellipsoids++;
        }
    }
}














// MultiEllipsoidSampler::drawFromHyperSphere()
//
// PURPOSE: 
//      Draws one point from the unit hyper-sphere and transforms
//      its coordinatesinto those of the input ellipsoid. The method is based on the 
//      approach described by Shaw J. R. et al. (2007; MNRAS, 378, 1365).
//
//
// INPUT:
//      ellipsoid: an object of class Ellipsoid that contains all the information
//      related to that ellipsoid.
//      drawnParametersPerObject: an Eigen Array of dimensions (Ndimensions) to contain the new 
//      coordinates for the point drawn from the unit sphere 
//      and transformed into coordinates of the input ellipsoid.
//
// OUTPUT:
//      void
//

void MultiEllipsoidSampler::drawFromHyperSphere(Ellipsoid &ellipsoid, RefArrayXd drawnParametersPerObject)
{
    // Pick a point uniformly from a unit hyper-sphere
    
    ArrayXd zeroCoordinates = ArrayXd::Zero(Ndimensions);
    double vectorNorm;

    do
    {
        for (int i = 0; i < Ndimensions; i++)
        {
            drawnParametersPerObject(i) = normal(engine);                            // Sample normally each coordinate
        }
        
        vectorNorm = metric.distance(drawnParametersPerObject,zeroCoordinates);
    }
    while (vectorNorm == 0);                                                // Repeat sampling if point falls in origin
        
    drawnParametersPerObject = drawnParametersPerObject/vectorNorm;                           // Normalize coordinates
    drawnParametersPerObject = pow(uniform(engine),1./Ndimensions)*drawnParametersPerObject;  // Sample uniformly in radial direction
    

    // Transform sphere coordinates to ellipsoid coordinates
    
    MatrixXd V = ellipsoid.getEigenvectorsMatrix().matrix();
    MatrixXd D = ellipsoid.getEigenvalues().sqrt().matrix().asDiagonal();
    MatrixXd T = V.transpose() * D;

    drawnParametersPerObject = (T * drawnParametersPerObject.matrix()) + ellipsoid.getCenterCoordinates().matrix();
}













// MultiEllipsoidSampler::intersection()
//
// PURPOSE:
//      Determines whether two input ellipsoids are overlapping according to
//      the algorithm described by Alfano & Greer (2003; Journal of Guidance,
//      Control and Dynamics, 26, 1).
//
// INPUT:
//      ellipsoid1: an object of class Ellipsoid that contains all the information
//      related to the first ellipsoid to be checked for intersection wtih the second ellispoid.
//      ellipsoid2: an object of class Ellipsoid that contains all the information
//      related to the second ellipsoid to be checked for intersection with the first ellipsoid.
//
// OUTPUT:
//      A boolean value specifying whether the two ellipsoids intersect (true) or not (false)
//
// REMARKS:
//      Coordinates of centers have to coincide in same order of the covariance matrix dimensions.
//      E.g. if first center coordinate is x, then first row and column in covariance matrix refer to x coordinate.
//

bool MultiEllipsoidSampler::intersection(Ellipsoid &ellipsoid1, Ellipsoid &ellipsoid2)
{
    // Construct translation matrix

    MatrixXd T1 = MatrixXd::Identity(Ndimensions+1,Ndimensions+1);
    MatrixXd T2 = MatrixXd::Identity(Ndimensions+1,Ndimensions+1);
    
    T1.bottomLeftCorner(1,Ndimensions) = -1.*ellipsoid1.getCenterCoordinates().transpose();
    T2.bottomLeftCorner(1,Ndimensions) = -1.*ellipsoid2.getCenterCoordinates().transpose();


    // Construct ellipsoid matrix in homogeneous coordinates

    MatrixXd A = MatrixXd::Zero(Ndimensions+1,Ndimensions+1);
    MatrixXd B = A;

    A(Ndimensions,Ndimensions) = -1;
    B(Ndimensions,Ndimensions) = -1;

    A.topLeftCorner(Ndimensions,Ndimensions) = ellipsoid1.getCovarianceMatrix().matrix().inverse();
    B.topLeftCorner(Ndimensions,Ndimensions) = ellipsoid2.getCovarianceMatrix().matrix().inverse();

    MatrixXd AT = T1*A*T1.transpose();        // Translating to ellispoid center
    MatrixXd BT = T2*B*T2.transpose();        // Translating to ellispoid center


    // Compute Hyper Quadric Matrix generating from the two ellipsoids 
    // and derive its eigenvalues decomposition

    MatrixXd C = AT.inverse() * BT;
    MatrixXcd CC(Ndimensions+1,Ndimensions+1);

    CC.imag() = MatrixXd::Zero(Ndimensions+1,Ndimensions+1); 
    CC.real() = C;
    
    ComplexEigenSolver<MatrixXcd> eigenSolver(CC);

    if (eigenSolver.info() != Success) abort();
    
    MatrixXcd E = eigenSolver.eigenvalues();
    MatrixXcd V = eigenSolver.eigenvectors();

    bool intersection = false;       // Start with no intersection
    double pointA;              // Point laying in elliposid A
    double pointB;              // Point laying in ellipsoid B
    
    for (int i = 0; i < Ndimensions+1; i++)      // Loop over all eigenvectors
    {
        if (V(Ndimensions,i).real() == 0)      // Skip inadmissible eigenvectors
            continue;                   
        else if (E(i).imag() != 0)
            {
                V.col(i) = V.col(i).array() * (V.conjugate())(Ndimensions,i);      // Multiply eigenvector by complex conjugate of last element
                V.col(i) = V.col(i).array() / V(Ndimensions,i).real();             // Normalize eigenvector to last component value
                pointA = V.col(i).transpose().real() * AT * V.col(i).real();        // Evaluate point from elliposid A
                pointB = V.col(i).transpose().real() * BT * V.col(i).real();        // Evaluate point from ellipsoid B

                if ((pointA <= 0) && (pointB <= 0))     // Accept only if point belongs to both ellipsoids
                {
                    intersection = true;                // Exit if intersection is found
                    break;
                }
            }
    }

    return intersection;

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
    overlappingEllipsoidsIndices.resize(1);
    nonOverlappingEllipsoidsIndices.resize(1);
    overlappingEllipsoidsIndices(0) = -1;           // Start with no overlapping ellipsoids found
    
    if (Nellipsoids == 1)         // If only one ellipsoid is found, then return one non-overlapping ellipsoid
    {
        nonOverlappingEllipsoidsIndices(0) = 0;
        return;
    }

    nonOverlappingEllipsoidsIndices(0) = -1;        // Start with no non-overlapping ellipsoids found

    bool overlapFlag = false;
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
            overlapFlag = intersection(ellipsoidsVector[i], ellipsoidsVector[j]);

            if (overlapFlag)    // If overlap occurred
            {
                countOverlap++;

                for (int k = 0; k < overlappingEllipsoidsIndices.size(); k++) // Check if j is saved already
                {
                    if (j == overlappingEllipsoidsIndices(k))       
                    {
                        saveFlagJ = false; // If j is saved already, don't save it again and go to next j
                        break;
                    }
                }

                
                if (saveFlagJ)      // If j is not saved yet, then save it
                {
                    overlappingEllipsoidsIndices.conservativeResize(countIndex+1);
                    overlappingEllipsoidsIndices(countIndex) = j;
                    countIndex++;
                }
            }
        }

        if (i == 0 && countOverlap != 0)   // If first i and at least one overlap is found, save also i and go to next i
        {
            overlappingEllipsoidsIndices.conservativeResize(countIndex+1);
            overlappingEllipsoidsIndices(countIndex) = i;
            countIndex++;
            continue;
        }
        else 
            if (i > 0 && countOverlap != 0)     // If i is not the first one and at least one overlap is found ...
            {
                for (int k = 0; k < overlappingEllipsoidsIndices.size(); k++) // Check if i is saved already
                {
                    if (i == overlappingEllipsoidsIndices(k))       
                    {
                        saveFlagI = false;       // If i is saved already, don't save it again and go to next i
                        break;
                    }
                }

                if (saveFlagI)      // If i is not saved yet, then save it
                {
                    overlappingEllipsoidsIndices.conservativeResize(countIndex+1);
                    overlappingEllipsoidsIndices(countIndex) = i;
                    countIndex++;
                }
            }
    }

    int Noverlaps;

    if (overlappingEllipsoidsIndices(0) != -1)
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
                if (i == overlappingEllipsoidsIndices(j))       // If ellipsoid i-th already overlap, then go to next i-th ellipsoid
                {
                    saveFlagI = false;
                    break;
                }
            }

            if (saveFlagI)      // If ellipsoid is not overlapping, save it among non-overlapping ellipsoids
            {   
                nonOverlappingEllipsoidsIndices(countIndex) = i;
                countIndex++;
            }
        }
    }

}



















// MultiEllipsoidSampler::pointIsInOverlap()
//
// PURPOSE:
//      Determines if the input point belongs to the input enlarged ellipsoid.
//
// INPUT:
//      ellipsoid: an Ellipsoid class object containing the information related
//      to the ellipsoid we want to check.
//      pointCoordinates: coordinates of the point to verify.
//
// OUTPUT:
//      A boolean value specifying whether the point belongs to the input
//      enlarged ellipsoid.
//

bool MultiEllipsoidSampler::pointIsInOverlap(Ellipsoid &ellipsoid, const RefArrayXd pointCoordinates)
{
    // Construct translation matrix

    MatrixXd T = MatrixXd::Identity(Ndimensions+1,Ndimensions+1);
    
    T.bottomLeftCorner(1,Ndimensions) = -1.*ellipsoid.getCenterCoordinates().transpose();


    // Construct ellipsoid matrix in homogeneous coordinates

    MatrixXd A = MatrixXd::Zero(Ndimensions+1,Ndimensions+1);
    A(Ndimensions,Ndimensions) = -1;
    
    MatrixXd V = ellipsoid.getEigenvectorsMatrix().matrix();
    MatrixXd C = MatrixXd::Zero(Ndimensions, Ndimensions);
    
    C = V * ellipsoid.getEigenvalues().matrix().asDiagonal() * V.transpose();      // Covariance matrix
    A.topLeftCorner(Ndimensions,Ndimensions) = C.inverse();

    MatrixXd AT = T*A*T.transpose();        // Translating to ellispoid center

    VectorXd X(Ndimensions+1);
    X.head(Ndimensions) = pointCoordinates.matrix();
    X(Ndimensions) = 1;

    bool overlap;
    overlap = false;        // Start with no overlap

    if (X.transpose() * AT * X <= 0)
        overlap = true;
        
    return overlap;
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

ArrayXi MultiEllipsoidSampler::getNonOverlappingEllipsoidsIndices()
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

ArrayXi MultiEllipsoidSampler::getOverlappingEllipsoidsIndices()
{
    return overlappingEllipsoidsIndices;
}














// MultiEllipsoidSampler::getEllipsoidsVector()
//
// PURPOSE:
//      Gets private data member ellipsoidsVector.
//
// OUTPUT:
//      a vector of Ellipsoids class objects containing the ellipsoids 
//      computed during the sampling process.
//

vector<Ellipsoid> MultiEllipsoidSampler::getEllipsoidsVector()
{
    return ellipsoidsVector;
}
