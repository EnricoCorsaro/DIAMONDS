#include "Ellipsoid.h"


// Ellipsoid::Ellipsoid()
//
// PURPOSE:
//      Class constructor.
//
// INPUT:
//      sample:                 Eigen Array of size (Ndimensions, sampleSize), 
//                              containing the coordinates of the points inside the ellipsoid.
//      enlargementFraction:    Initial enlargement fraction to compute the ellipsoids for the first time
//

Ellipsoid::Ellipsoid(RefArrayXXd sample, const double enlargementFraction)
: sample(sample),
  sampleSize(sample.cols()),
  Ndimensions(sample.rows()),
  uniform(0.0, 1.0),
  normal(0.0, 1.0)
{
    // Set the seed of the random generator using the clock

    clock_t clockticks = clock();
    engine.seed(clockticks);


    // Resize the matrices to their proper size

    originalEigenvalues.resize(Ndimensions);
    enlargedEigenvalues.resize(Ndimensions);
    centerCoordinates.resize(Ndimensions);
    eigenvectors.resize(Ndimensions, Ndimensions);
    covarianceMatrix.resize(Ndimensions, Ndimensions);


    // Compute the covariance matrix of the sample of points

    Functions::clusterCovariance(sample, covarianceMatrix, centerCoordinates);
   
    
    // Compute the eigenvectors and eigenvalues of this covariance matrix

    bool covarianceDecompositionIsSuccessful = Functions::selfAdjointMatrixDecomposition(covarianceMatrix, originalEigenvalues, eigenvectors);
    if (!covarianceDecompositionIsSuccessful) abort();


    // Set the enlargement fraction, and (re)compute the corresponding eigenvalues, eigenvectors etc.

    resetEnlargementFraction(enlargementFraction);
}











// Ellipsoid::~Ellipsoid()
//
// PURPOSE:
//      Class destructor.
//

Ellipsoid::~Ellipsoid()
{

}











// Ellipsoid::resetEnlargementFraction()
//
// PURPOSE: 
//      Compute covariance matrix, center coordinates, 
//      eigenvalues and eigenvectors matrix of the ellipsoid.
//      Eigenvalues are automatically enlarged according to the
//      enlargementFraction term.
//
// INPUT:
//      newEnlargementFraction:  a double to contain the new enlargement
//                               fraction to be used for the chosen ellipsoid.
//
// OUTPUT:
//      void

void Ellipsoid::resetEnlargementFraction(const double newEnlargementFraction)
{
    // Save the new user-specified enlargement fraction

    this->enlargementFraction = newEnlargementFraction;
   

    // Enlarge the eigenvalues with the user-specified fraction

    enlargedEigenvalues = originalEigenvalues.sqrt() + enlargementFraction * originalEigenvalues.sqrt();
    enlargedEigenvalues = enlargedEigenvalues.square();


    // Recompute the hypervolume contained in the enlarged ellipsoid

    hyperVolume = enlargedEigenvalues.sqrt().prod();


    // Recompute the covariance matrix with the enlarged eigenvalues

    covarianceMatrix = eigenvectors.matrix() * enlargedEigenvalues.matrix().asDiagonal() * eigenvectors.matrix().transpose();
}











// Ellipsoid::overlapsWith()
//
// PURPOSE:
//      Determines whether this ellipsoid overlaps with the given ellipsoid.
//      The algorithm used is the one described by Alfano & Greer (2003; 
//      Journal of Guidance, Control and Dynamics, 26, 1).
//
// INPUT:
//      ellipsoid:                                  Ellipsoid object
//      ellipsoidMatrixDecompositionIsSuccessful:   a boolean specifying if the decomposition of 
//                                                  the ellipsoid matrix was successful
//
// OUTPUT:
//      A boolean value specifying whether the two ellipsoids overlap (true) or not (false)
//
// REMARKS:
//      The coordinate system used by the ellipsoid in the argument must be the same as the
//      one used by this ellipsoid. This is relevant to know which coordinate corresponds to
//      which column in the covariance matrices of the ellipsoids.
//

bool Ellipsoid::overlapsWith(Ellipsoid ellipsoid, bool &ellipsoidMatrixDecompositionIsSuccessful)
{
    // Construct translation matrix

    MatrixXd T1 = MatrixXd::Identity(Ndimensions+1,Ndimensions+1);
    MatrixXd T2 = MatrixXd::Identity(Ndimensions+1,Ndimensions+1);
    
    T1.bottomLeftCorner(1,Ndimensions) = (-1.0) * centerCoordinates.transpose();
    T2.bottomLeftCorner(1,Ndimensions) = (-1.0) * ellipsoid.getCenterCoordinates().transpose();


    // Construct ellipsoid matrix in homogeneous coordinates

    MatrixXd A = MatrixXd::Zero(Ndimensions+1,Ndimensions+1);
    MatrixXd B = A;

    A(Ndimensions,Ndimensions) = -1;
    B(Ndimensions,Ndimensions) = -1;

    A.topLeftCorner(Ndimensions,Ndimensions) = covarianceMatrix.matrix().inverse();
    B.topLeftCorner(Ndimensions,Ndimensions) = ellipsoid.getCovarianceMatrix().matrix().inverse();

    MatrixXd AT = T1*A*T1.transpose();        // Translating to ellipsoid center
    MatrixXd BT = T2*B*T2.transpose();        // Translating to ellipsoid center


    // Compute Hyper Quadric Matrix generating from the two ellipsoids 
    // and derive its eigenvalues decomposition

    MatrixXd C = AT.inverse() * BT;
    MatrixXcd CC(Ndimensions+1,Ndimensions+1);

    CC.imag() = MatrixXd::Zero(Ndimensions+1,Ndimensions+1); 
    CC.real() = C;
    
    ComplexEigenSolver<MatrixXcd> eigenSolver(CC);


    // If eigenvalue decomposition fails, set control flag to false 
    // to stop the nested sampling and print the results 

    if (eigenSolver.info() != Success)
    {
        ellipsoidMatrixDecompositionIsSuccessful = false;
    }
    
    MatrixXcd E = eigenSolver.eigenvalues();
    MatrixXcd V = eigenSolver.eigenvectors();

    bool ellipsoidsDoOverlap = false;       // Assume no overlap in the beginning
    double pointA;                          // Point laying in this ellipsoid
    double pointB;                          // Point laying in the other ellipsoid


    // Loop over all eigenvectors

    for (int i = 0; i < Ndimensions+1; i++) 
    {
        // Skip inadmissible eigenvectors

        if (V(Ndimensions,i).real() == 0)
        {
            continue;                   
        }
        else if (E(i).imag() != 0)
            {
                V.col(i) = V.col(i).array() * (V.conjugate())(Ndimensions,i);      // Multiply eigenvector by complex conjugate of last element
                V.col(i) = V.col(i).array() / V(Ndimensions,i).real();             // Normalize eigenvector to last component value
                pointA = V.col(i).transpose().real() * AT * V.col(i).real();       // Evaluate point from this ellipsoid
                pointB = V.col(i).transpose().real() * BT * V.col(i).real();       // Evaluate point from the other ellipsoid


                // Accept only if point belongs to both ellipsoids

                if ((pointA <= 0) && (pointB <= 0))  
                {
                    ellipsoidsDoOverlap = true;            // Exit if ellipsoidsDoOverlap is found
                    break;
                }
            }
    }

    return ellipsoidsDoOverlap;
}










// Ellipsoid::containsPoint()
//
// PURPOSE:
//      Determines if this ellipsoid contains the given point
//
// INPUT:
//      pointCoordinates: coordinates of the point to verify.
//
// OUTPUT:
//      'true' if the point falls within the (possibly enlarged) boundaries
//      of this ellipsoid, 'false' otherwise.
//

bool Ellipsoid::containsPoint(const RefArrayXd pointCoordinates)
{
    // Construct translation matrix

    MatrixXd T = MatrixXd::Identity(Ndimensions+1,Ndimensions+1);
    
    T.bottomLeftCorner(1,Ndimensions) = (-1.) * centerCoordinates.transpose();


    // Construct ellipsoid matrix in homogeneous coordinates

    MatrixXd A = MatrixXd::Zero(Ndimensions+1,Ndimensions+1);
    A(Ndimensions,Ndimensions) = -1;
    
    MatrixXd C = MatrixXd::Zero(Ndimensions, Ndimensions);


    // Compute the covariance matrix

    C =  eigenvectors.matrix() * enlargedEigenvalues.matrix().asDiagonal() 
                               * eigenvectors.matrix().transpose(); 
    A.topLeftCorner(Ndimensions,Ndimensions) = C.inverse();


    // Translate to the ellipsoid center

    MatrixXd AT = T * A * T.transpose(); 

    VectorXd X(Ndimensions+1);
    X.head(Ndimensions) = pointCoordinates.matrix();
    X(Ndimensions) = 1;


    // Check if the point belongs to this ellipsoid

    bool pointBelongsToThisEllipsoid;

    if (X.transpose() * AT * X <= 0)
    {
        pointBelongsToThisEllipsoid = true;
    }
    else
    {
        pointBelongsToThisEllipsoid = false;
    }
        
    return pointBelongsToThisEllipsoid;
}









// Ellipsoid::drawPoint()
//
// PURPOSE: 
//      Draw a random point inside this ellipsoid. The algorithm used is the one 
//      described by Shaw J. R. et al. (2007; MNRAS, 378, 1365). First a point is
//      drawn from the unit hypersphere, whose coordinates are then transformed
//      to land into the ellipsoid.
//
//
// INPUT:
//      drawnPoint: Eigen Array that will contain the N-dimensional coordinates of the
//                  newly drawn point.
//
// OUTPUT:
//      void
//

void Ellipsoid::drawPoint(RefArrayXd drawnPoint)
{
    // Pick a point uniformly from a unit hyper-sphere
    
    do
    {
        for (int i = 0; i < Ndimensions; i++)
        {
            // Sample normally in each coordinate direction

            drawnPoint(i) = normal(engine); 
        }
    }
    while ((drawnPoint == 0.0).all());    // Repeat sampling if point falls in origin
    

    // Normalize the point so that it belongs to the unit hyper-sphere
    
    drawnPoint = drawnPoint / drawnPoint.matrix().norm(); 
    
 

    // Sample uniformly in the radial direction

    drawnPoint = pow(uniform(engine), 1./Ndimensions) * drawnPoint; 
    

    // Transform sphere coordinates to ellipsoid coordinates
    
    MatrixXd D = enlargedEigenvalues.sqrt().matrix().asDiagonal();
    MatrixXd T = eigenvectors.matrix() * D;    
    
    drawnPoint = (T * drawnPoint.matrix()) + centerCoordinates.matrix();
}










// Ellipsoid::getCenterCoordinates()
//
// PURPOSE: 
//      Gets the protected data member centerCoordinates.      
//
// OUTPUT:
//      An Eigen Array of dimensions (Ndimensions)
//      containing all the center coordinates of the ellipsoid.
//

ArrayXd Ellipsoid::getCenterCoordinates()
{
    return centerCoordinates;
}











// Ellipsoid::getEigenvalues()
//
// PURPOSE: 
//      Gets the protected data member eigenvalues.      
//
// OUTPUT:
//      An Eigen Array of dimensions (Ndimensions), containing all enlarged eigenvalues of the ellipsoid.
//

ArrayXd Ellipsoid::getEigenvalues()
{
    return enlargedEigenvalues;
}












// Ellipsoid::getSample()
//
// PURPOSE: 
//      Gets the protected data member sample.      
//
// OUTPUT:
//      An Eigen Array matrix of dimensions (Ndimensions, sampleSize) 
//      containing all the coordinates of the objects contained in the
//      ellipsoid.
//

ArrayXXd Ellipsoid::getSample()
{
    return sample;
}












// Ellipsoid::getCovarianceMatrix()
//
// PURPOSE: 
//      Gets the protected data member enlargedCovarianceMatrix.      
//
// OUTPUT:
//      An Eigen Array matrix of dimensions (Ndimensions, Ndimensions) 
//      containing the (enlarged) covariance matrix of the ellipsoid.
//

ArrayXXd Ellipsoid::getCovarianceMatrix()
{
    return covarianceMatrix;
}












// Ellipsoid::getEigenvectors()
//
// PURPOSE: 
//      Gets the protected data member enlargedEigenvectors.      
//
// OUTPUT:
//      An Eigen Array matrix of dimensions (Ndimensions, Ndimensions), 
//      containing all eigenvectors of the ellipsoid.
//

ArrayXXd Ellipsoid::getEigenvectors()
{
    return eigenvectors;
}












// Ellipsoid::getSampleSize()
//
// PURPOSE: 
//      Gets the number of points of the sample inside the ellipsoid
//
// OUTPUT:
//      An integer 
//

int Ellipsoid::getSampleSize()
{
    return sampleSize;
}












// Ellipsoid::getHyperVolume()
//
// PURPOSE: 
//      Get the private data member hyperVolumes.      
//
// OUTPUT:
//      A double containing the hyper-volume of the ellipsoid.
//

double Ellipsoid::getHyperVolume()
{
    return hyperVolume;
}












// Ellipsoid::getEnlargementFraction()
//
// PURPOSE: 
//      Get the private data member enlargementFraction.      
//
// OUTPUT:
//      A double containing the enlargement fraction computed for the ellipsoid
//      and updated to the last resetting.
//

double Ellipsoid::getEnlargementFraction()
{
    return enlargementFraction;
}
