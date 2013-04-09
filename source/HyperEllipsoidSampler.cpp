#include "HyperEllipsoidSampler.h"

// HyperEllipsoidSampler::HyperEllipsoidSampler()
//
// PURPOSE: 
//      Base class constructor.
//
// INPUT:
//      prior: prior object containing the prior distributions for each parameter
//      likelihood: likelihood object containing the likelihood and the model to be used
//      metric: metric object to define the distance between two points
//      Nobjects: the initial number of objects coming from the main nested sampling loop
//      initialEnlargementFactor: initial value of the enlargement for the ellipsoids
//      alpha: a value between 0 and 1 that defines the shrinking rate of the enlargement factors
//

HyperEllipsoidSampler::HyperEllipsoidSampler(Prior &prior, Likelihood &likelihood, Metric &metric, const int Nobjects, const double initialEnlargementFactor, const double alpha)
: HyperQuadricSampler(prior, likelihood, metric, Nobjects),
  engine(time(0)),
  uniform(0.0,1.0),
  normal(0.0,1.0),
  initialEnlargementFactor(initialEnlargementFactor),
  alpha(alpha)
{

}











// HyperEllipsoidSampler::~HyperEllipsoidSampler()
//
// PURPOSE: 
//      Base class destructor.
//

HyperEllipsoidSampler::~HyperEllipsoidSampler()
{

}












// HyperEllipsoidSampler::drawWithConstraint()
//
// PURPOSE:
//      Draws a new object from the largest ellipsoid found in the sample.
//
// INPUT:
//      totalSampleOfParameters: an Eigen Array matrix of size (Ndimensions, Nobjects)
//      Nclusters: optimal number of clusters found by clustering algorithm
//      clusterIndices: indices of clusters for each point of the sample
//      logWidthInPriorMass: log Value of prior volume from the actual nested iteration.
//      nestedSampleOfParameters: an Eigen Array matrix of size (Ndimensions,Ndraws) to contain the
//      coordinates of the drawn point to be used for the next nesting loop. When used
//      for the first time, the array contains the coordinates of the worst nested
//      object to be updated.
//      likelihood: an object of class Likelihood to contain the likelihood function
//      used in the sampling criterion.
//      
//
// OUTPUT:
//      void
//

void HyperEllipsoidSampler::drawWithConstraint(const RefArrayXXd totalSampleOfParameters, const int Nclusters, const RefArrayXi clusterIndices, const double logWidthInPriorMass, RefArrayXXd nestedSampleOfParameters)
{
    assert(totalSampleOfParameters.cols() == clusterIndices.size());
    assert(nestedSampleOfParameters.rows() == totalSampleOfParameters.rows());
    assert(Nclusters > 0);

    int Ndimensions = totalSampleOfParameters.rows();
    int Ndraws = nestedSampleOfParameters.cols();
    HyperEllipsoidIntersector intersector;


    // Compute covariance matrices and centers for ellipsoids associated with each cluster.
    // Compute eigenvalues and eigenvectors, enlarge eigenvalues and compute their hyper volumes.

    computeEllipsoids(totalSampleOfParameters, Nclusters, clusterIndices, logWidthInPriorMass);
    
    
    // Find which ellipsoids are overlapping and which are not
    
    intersector.findOverlappingEllipsoids(Nclusters, allEnlargedEigenvalues, allEigenvectorsMatrix, allCentersCoordinates);
    ArrayXi nonOverlappingEllipsoidsIndices = intersector.getNonOverlappingEllipsoidsIndices();
    ArrayXi overlappingEllipsoidsIndices = intersector.getOverlappingEllipsoidsIndices();
 
    
    // Isolated ellipsoid sampling

    if (nonOverlappingEllipsoidsIndices(0) != -1)
    {
        // If non-overlapping ellipsoids exist start separate nested sampling 
        // from each non overlapping enlarged ellipsoid 
    
        int clusterIndex;       // Cluster index
        int NobjectsInCluster;          // Number of objects for the nested sampling within the identified ellipsoid
        double logLocalWidthInPriorMass;    // Shrinked prior volume for separate nesting of ellipsoid

        for (int n = 0; n < nonOverlappingEllipsoidsIndices.size(); n++)
        {
            clusterIndex = nonOverlappingEllipsoidsIndices(n);
            NobjectsInCluster = NobjectsPerCluster(clusterIndex);
            logLocalWidthInPriorMass = exp(logWidthInPriorMass)*(NobjectsInCluster/Nobjects);
            cout << "Start separate nesting loop for ellipsoid #" << clusterIndex << endl;
            // NestedSampler nestedSampler(prior, likelihood);          // TO BE FIXED
            // nestedSampler.runSubordinate(NobjectsInCluster,clusterSampleOfObjects,logLocalWidthInPriorMass);
        }
    }


    // Overlapping ellipsoid sampling

    if (overlappingEllipsoidsIndices(0) != -1)
    {
        // If overlapping ellipsoids exist choose one with probability according to its volume fraction

        int NoverlappingEllipsoids = overlappingEllipsoidsIndices.size();
        int ellipsoidIndex;
        double volumeProbability;
        double actualProbability; 
        double totalVolume;
        uniform_int_distribution<int> uniformIndex(0,NoverlappingEllipsoids-1);

        for (int m = 0; m < NoverlappingEllipsoids; m++)
        {   
            ellipsoidIndex = overlappingEllipsoidsIndices(m);
            totalVolume += hyperVolumes(ellipsoidIndex);                   // Compute total volume of overlapping ellipsoids
        }

        do
        {
            ellipsoidIndex = overlappingEllipsoidsIndices(uniformIndex(engine));      // Pick up one ellipsoid randomly
            volumeProbability = hyperVolumes(ellipsoidIndex)/totalVolume;             // Compute the probability for the selected ellipsoid based on its volume
            actualProbability = uniform(engine);            // Give a probability value between 0 and 1
        }
        while (actualProbability > volumeProbability);      // If actualProbability < volumeProbability then pick the corresponding ellipsoid

        ArrayXd centerCoordinates1;
        ArrayXd ellipsoidEigenvalues1;
        ArrayXXd ellipsoidEigenvectorsMatrix1;
        centerCoordinates1 = allCentersCoordinates.segment(ellipsoidIndex*Ndimensions, Ndimensions);
        ellipsoidEigenvalues1 = allEnlargedEigenvalues.segment(ellipsoidIndex*Ndimensions, Ndimensions);
        ellipsoidEigenvectorsMatrix1 = allEigenvectorsMatrix.block(0, ellipsoidIndex*Ndimensions, Ndimensions, Ndimensions);


        // Create arrays for a second ellipsoid to check for overlap

        int ellipsoidIndex2;
        double rejectProbability;
        double logLikelihood;
        ArrayXd centerCoordinates2 = ArrayXd::Zero(Ndimensions);
        ArrayXd ellipsoidEigenvalues2 = ArrayXd::Zero(Ndimensions);
        ArrayXXd ellipsoidEigenvectorsMatrix2 = ArrayXXd::Zero(Ndimensions, Ndimensions);
        ArrayXd drawnParameters;

        
        // Sampe uniformly from first chosen ellipsoid until new parameter with logLikelihood > logLikelihoodConstraint is found
        // If parameter is contained in Noverlaps ellipsoids, then accept parameter with probability 1/Noverlaps

        //double logLikelihoodConstraint = likelihood.logValue(drawnParameters);
        
        for (int m = 0; m < Ndraws; m++)
        {
            drawnParameters = nestedSampleOfParameters.col(m);
            
            do
            {
                //do
                //{
                    drawFromHyperSphere(ellipsoidEigenvalues1, ellipsoidEigenvectorsMatrix1, centerCoordinates1, drawnParameters);
                //   logLikelihood = likelihood.logValue(drawnParameters);
                //}
                //while (logLikelihood <= logLikelihoodConstraint);

                int Noverlaps = 1;        // Start with self-overlap only
                bool overlapFlag = false;

                for (int i = 0; i < NoverlappingEllipsoids; i++)        // Check other possibile overlapping ellipsoids
                {
                    if (overlappingEllipsoidsIndices(i) == ellipsoidIndex)     // Skip if self-overlap
                        continue;
                    else
                    {
                        ellipsoidIndex2 = overlappingEllipsoidsIndices(i);
                        centerCoordinates2 = allCentersCoordinates.segment(ellipsoidIndex2*Ndimensions, Ndimensions);
                        ellipsoidEigenvalues2 = allEnlargedEigenvalues.segment(ellipsoidIndex2*Ndimensions, Ndimensions);
                        ellipsoidEigenvectorsMatrix2 = allEigenvectorsMatrix.block(0, ellipsoidIndex2*Ndimensions, Ndimensions, Ndimensions);


                        // Check if point belongs to ellipsoid 2

                        overlapFlag = intersector.checkPointForOverlap(ellipsoidEigenvalues2, ellipsoidEigenvectorsMatrix2, centerCoordinates2, drawnParameters);
                
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
       
            nestedSampleOfParameters.col(m) = drawnParameters;  // Save set of parameters for one object into total sample
        }
    } // END sampling from overlapping ellipsoids
}












// HyperEllipsoidSampler::computeEllipsoids()
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
//      logWidthInPriorMass: log Value of prior volume from the actual nested iteration.
//
// OUTPUT:
//      void
//

void HyperEllipsoidSampler::computeEllipsoids(const RefArrayXXd totalSampleOfParameters, const int Nclusters, const RefArrayXi clusterIndices, const double logWidthInPriorMass)
{
    assert(totalSampleOfParameters.cols() == clusterIndices.size());
    int Ndimensions = totalSampleOfParameters.rows();       

    allClustersCovarianceMatrix = ArrayXXd::Zero(Ndimensions, Nclusters*Ndimensions);
    allCentersCoordinates = ArrayXd::Zero(Nclusters * Ndimensions);
    NobjectsPerCluster = ArrayXi::Zero(Nclusters);
    
    ArrayXXd covarianceMatrix(Ndimensions, Ndimensions);
    ArrayXd centerCoordinates(Ndimensions);
    
    ArrayXd ellipsoidEigenvalues(Ndimensions);
    ArrayXXd ellipsoidEigenvectorsMatrix(Ndimensions, Ndimensions);
    
    hyperVolumes.resize(Nclusters);                                            // Array containing hyper-volumes of each ellipsoid
    allEigenvectorsMatrix.resize(Ndimensions, Nclusters*Ndimensions);          // Matrix containing all eigenvectors of each ellipsoid
    allEnlargedEigenvalues.resize(Nclusters*Ndimensions);                      // Array containing all enlarged eigenvalues of each ellipsoid 
    allEigenvalues.resize(Nclusters*Ndimensions);                              // Array containing all original eigenvalues of each ellipsoid 


    // Divide the sample according to the clustering done

    for (int i = 0; i < Nobjects; i++)         // Find number of points (objects) per cluster
    {
        NobjectsPerCluster(clusterIndices(i))++;
    }

    ArrayXd oneDimensionSampleOfParameters(Nobjects);
    ArrayXXd totalSampleOfParametersOrdered(Ndimensions, Nobjects);
    ArrayXi clusterIndicesCopy(Nobjects);
    clusterIndicesCopy = clusterIndices;

    for (int i = 0; i < Ndimensions; i++)     // Order points in each dimension according to cluster indices
    {
        oneDimensionSampleOfParameters = totalSampleOfParameters.row(i);
        Functions::sortElementsInt(clusterIndicesCopy, oneDimensionSampleOfParameters);
        totalSampleOfParametersOrdered.row(i) = oneDimensionSampleOfParameters;
        clusterIndicesCopy = clusterIndices;
    }


    // Compute covariance matrices and center coordinates for different blocks corresponding to different clusters
    // Then perform eigenvalues decomposition, enlarge eigenvalues and save enlarged eigenvalues and eigenvectors
    
    ArrayXXd clusterSample;
    int actualNobjects = 0;

    for (int i = 0; i < Nclusters; i++)
    {   
        if (NobjectsPerCluster(i) <= Ndimensions + 1)        // Skip cluster if number of points is not large enough
        {
            actualNobjects += NobjectsPerCluster(i);
            continue;
        }

        clusterSample = ArrayXXd::Zero(Ndimensions, NobjectsPerCluster(i));
        clusterSample = totalSampleOfParametersOrdered.block(0, actualNobjects, Ndimensions, NobjectsPerCluster(i));
        actualNobjects += NobjectsPerCluster(i);
        
        Functions::clusterCovariance(clusterSample, covarianceMatrix, centerCoordinates);
        
        allClustersCovarianceMatrix.block(0, i*Ndimensions, Ndimensions, Ndimensions) = covarianceMatrix; 
        allCentersCoordinates.segment(i*Ndimensions, Ndimensions) = centerCoordinates;
        
        Functions::selfAdjointMatrixDecomposition(covarianceMatrix, ellipsoidEigenvalues, ellipsoidEigenvectorsMatrix); // Do eigenvalue decomposition
        allEigenvalues.segment(i*Ndimensions, Ndimensions) = ellipsoidEigenvalues;
        ellipsoidEnlarger(ellipsoidEigenvalues, logWidthInPriorMass, NobjectsPerCluster(i));                        // Enlarge eigenvalues
        
        hyperVolumes(i) = sqrt(ellipsoidEigenvalues.prod());                                                       // Compute hyper-volume
        allEigenvectorsMatrix.block(0, i*Ndimensions, Ndimensions, Ndimensions) = ellipsoidEigenvectorsMatrix;
        allEnlargedEigenvalues.segment(i*Ndimensions, Ndimensions) = ellipsoidEigenvalues;
    }

}












// HyperEllipsoidSampler::getAllClustersCovarianceMatrix()
//
// PURPOSE: 
//      Gets the private data member allClustersCovarianceMatrix.      
//
// OUTPUT:
//      An Eigen Array matrix of dimensions (Ndimensions, Nclusters*Ndimensions) 
//      partitioned in blocks of dimensions (Ndimensions, Ndimensions) each, 
//      containing all the covariance matrices of the clusters identified.
//

ArrayXXd HyperEllipsoidSampler::getAllClustersCovarianceMatrix()
{
    return allClustersCovarianceMatrix;
}













// HyperEllipsoidSampler::getAllCentersCoordinates()
//
// PURPOSE: 
//      Gets the private data member allCentersCoordinates.      
//
// OUTPUT:
//      An Eigen Array of dimensions (Nclusters*Ndimensions)
//      containing all the arrays of the center coordinates of each cluster.
//

ArrayXd HyperEllipsoidSampler::getAllCentersCoordinates()
{
    return allCentersCoordinates;
}












// HyperEllipsoidSampler::getNobjectsPerCluster()
//
// PURPOSE: 
//      Gets the private data member NobjectsPerCluster.      
//
// OUTPUT:
//      An Eigen Array of integers of dimensions (Nclusters)
//      containing the number of points in each cluster.
//

ArrayXi HyperEllipsoidSampler::getNobjectsPerCluster()
{
    return NobjectsPerCluster;
}













// HyperEllipsoidSampler::getAllEigenvectorsMatrix()
//
// PURPOSE: 
//      Gets the private data member allEnlargedEigenvectors.      
//
// OUTPUT:
//      An Eigen Array matrix of dimensions 
//      (Ndimensions, Nclusters * Ndimensions), containing 
//      all eigenvectors of each ellipsoid.
//

ArrayXXd HyperEllipsoidSampler::getAllEigenvectorsMatrix()
{
    return allEigenvectorsMatrix;
}














// HyperEllipsoidSampler::getAllEigenvalues()
//
// PURPOSE: 
//      Gets the private data member allEigenvalues.      
//
// OUTPUT:
//      An Eigen Array of dimensions (Nclusters * Ndimensions),
//      containing all original eigenvalues of each ellipsoid.
//

ArrayXd HyperEllipsoidSampler::getAllEigenvalues()
{
    return allEigenvalues;
}














// HyperEllipsoidSampler::getAllEnlargedEigenvalues()
//
// PURPOSE: 
//      Gets the private data member allEnlargedEigenvalues.      
//
// OUTPUT:
//      An Eigen Array of dimensions (Nclusters * Ndimensions),
//      containing all enlarged eigenvalues of each ellipsoid.
//

ArrayXd HyperEllipsoidSampler::getAllEnlargedEigenvalues()
{
    return allEnlargedEigenvalues;
}












// HyperEllipsoidSampler::getAllEnlargedCovarianceMatrix()
//
// PURPOSE: 
//      Computes the private data member allEnlargedCovarianceMatrix and returns it
//      as an output.      
//
// OUTPUT:
//      An Eigen Array matrix of dimensions (Ndimensions, Nclusters*Ndimensions) 
//      partitioned in blocks of dimensions (Ndimensions, Ndimensions) each, 
//      containing all the covariance matrices of the clusters identified
//      computed with enlarged eigenvalues.
//

ArrayXXd HyperEllipsoidSampler::getAllEnlargedCovarianceMatrix()
{
    computeAllEnlargedCovarianceMatrix();

    return allEnlargedCovarianceMatrix;
}
















// HyperEllipsoidSampler::getHyperVolumes()
//
// PURPOSE: 
//      Gets the private data member hyperVolumes.      
//
// OUTPUT:
//      An Eigen Array of dimensions (Nclusters),
//      containing the hyper-volume each ellipsoid.
//

ArrayXd HyperEllipsoidSampler::getHyperVolumes()
{
    return hyperVolumes;
}












// HyperEllipsoidSampler::ellipsoidEnlarger()
//
// PURPOSE: 
//      Enlarge semi-major axes of the input ellipsoid according to the
//      recepie given in Feroz & Hobson (2008; MNRAS, 384, 449).
//
// INPUT:
//      ellipsoidEigenvalues: an Eigen Array containing the eigenvalues to be enlarged. 
//      logWidthInPriorMass: log Value of prior volume from the actual nested iteration.
//      NobjectsInCluster: number of points used to compute the ellipsoid.
//      of the associated ellipsoid covariance matrix.
//
// OUTPUT:
//      void
//

void HyperEllipsoidSampler::ellipsoidEnlarger(RefArrayXd ellipsoidEigenvalues, const double logWidthInPriorMass, const int NobjectsInCluster)
{
    double enlargementFactor;
    enlargementFactor = initialEnlargementFactor * pow(exp(logWidthInPriorMass), alpha) * sqrt(Nobjects/NobjectsInCluster);
  

    // Enlarge semi-major axes of the ellispoid
    
    ArrayXd D = ellipsoidEigenvalues.sqrt() + enlargementFactor*ellipsoidEigenvalues.sqrt();
   

    // Adjust eigenvalues with enlarged semi-major axes

    ellipsoidEigenvalues = D * D;

}











// HyperEllipsoidSampler::drawFromHyperSphere()
//
// PURPOSE: 
//      Draws one point from the unit hyper-sphere and transforms
//      its coordinatesinto those of the input ellipsoid. The method is based on the 
//      approach described by Shaw J. R. et al. (2007; MNRAS, 378, 1365).
//
//
// INPUT:
//      ellipsoidEigenvalues: an Eigen Array containing the eigenvalues 
//      of the associated ellipsoid covariance matrix.
//      eigenVectorsMatrix: an Eigen Array containing the eigen vectors of the 
//      covariance matrix as column vectors.
//      centerCoordinates: an Eigen Array containing the coordinates of the 
//      central point of the ellipsoid.
//      drawnParameters: an Eigen Array of dimensions (Ndimensions) to contain the new 
//      coordinates for the point drawn from the unit sphere 
//      and transformed into coordinates of the input ellipsoid.
//
// OUTPUT:
//      void
//

void HyperEllipsoidSampler::drawFromHyperSphere(const RefArrayXd ellipsoidEigenvalues, const RefArrayXXd eigenVectorsMatrix, const RefArrayXd centerCoordinates, RefArrayXd drawnParameters)
{
    assert(eigenVectorsMatrix.rows() == eigenVectorsMatrix.cols());
    assert(ellipsoidEigenvalues.size() == centerCoordinates.size());
    assert(drawnParameters.size() == ellipsoidEigenvalues.size());
    assert(eigenVectorsMatrix.rows() == centerCoordinates.size());


    // Pick a point uniformly from a unit hyper-sphere
    
    int Ndimensions = centerCoordinates.size();
    ArrayXd zeroCoordinates = ArrayXd::Zero(Ndimensions);
    double vectorNorm;

    do
    {
        for (int i = 0; i < Ndimensions; i++)
        {
            drawnParameters(i) = normal(engine);            // Sample normally each coordinate
        }
        
        vectorNorm = metric.distance(drawnParameters,zeroCoordinates);
    }
    while (vectorNorm == 0);                            // Repeat sampling if point falls in origin
        
    drawnParameters = drawnParameters/vectorNorm;       // Normalize coordinates
    drawnParameters = pow(uniform(engine),1./Ndimensions)*drawnParameters;  // Sample uniformly in radial direction
    

    // Transform sphere coordinates to ellipsoid coordinates

    MatrixXd V = eigenVectorsMatrix.matrix();
    MatrixXd D = ellipsoidEigenvalues.sqrt().matrix().asDiagonal();
    MatrixXd T = V.transpose() * D;

    drawnParameters = (T * drawnParameters.matrix()) + centerCoordinates.matrix();
}













// HyperEllipsoidSampler::computeAllEnlargedCovarianceMatrix()
// 
// PURPOSE:
//      Compute the total covariance matrix corresponding to the
//      enlarged eigenvalues of each ellipsoid.
//  
// OUTPUT:
//      void
//

void HyperEllipsoidSampler::computeAllEnlargedCovarianceMatrix()
{
    int Ndimensions = allEigenvectorsMatrix.rows();
    int Nclusters = allEnlargedEigenvalues.size()/Ndimensions;
    ArrayXd ellipsoidEigenvalues(Ndimensions);
    ArrayXXd ellipsoidEigenvectorsMatrix(Ndimensions,Ndimensions);
    ArrayXXd covarianceMatrix(Ndimensions,Ndimensions);
    allEnlargedCovarianceMatrix.resize(Ndimensions,Ndimensions*Nclusters);

    for (int i=0; i < Nclusters; i++)
    {
        ellipsoidEigenvalues = allEnlargedEigenvalues.segment(i*Ndimensions, Ndimensions);
        ellipsoidEigenvectorsMatrix = allEigenvectorsMatrix.block(0, i*Ndimensions, Ndimensions, Ndimensions);
        covarianceMatrix = ellipsoidEigenvectorsMatrix.matrix() * ellipsoidEigenvalues.matrix().asDiagonal() * ellipsoidEigenvectorsMatrix.matrix().transpose();
        allEnlargedCovarianceMatrix.block(0, i*Ndimensions, Ndimensions, Ndimensions) = covarianceMatrix;
    }
}
