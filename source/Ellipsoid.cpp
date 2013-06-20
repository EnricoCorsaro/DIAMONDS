#include "Ellipsoid.h"


// Ellipsoid::Ellipsoid()
//
// PURPOSE:
//      Class constructor.
//
// INPUT:
//      sampleOfParameters: an Eigen Array matrix of size (Ndimensions, Nobjects), 
//                          containing the coordinates of the objects inside the ellipsoid.
//      enlargementFactor: TODO
//

Ellipsoid::Ellipsoid(RefArrayXXd sampleOfParameters, const double enlargementFactor)
: sampleOfParameters(sampleOfParameters),
  Nobjects(sampleOfParameters.cols()),
  Ndimensions(sampleOfParameters.rows())
{
    // Resize the matrices to their proper size

    originalEigenvalues.resize(Ndimensions);
    enlargedEigenvalues.resize(Ndimensions);
    centerCoordinates.resize(Ndimensions);
    eigenvectors.resize(Ndimensions, Ndimensions);
    covarianceMatrix.resize(Ndimensions, Ndimensions);

    // Compute the covariance matrix of the sample of points

    Functions::clusterCovariance(sampleOfParameters, covarianceMatrix, centerCoordinates);

    // Compute the eigenvectors and eigenvalues of this covariance matrix

    Functions::selfAdjointMatrixDecomposition(covarianceMatrix, originalEigenvalues, eigenvectors);

    // Set the enlargement factor, and (re)compute the corresponding eigenvalues, eigenvectors etc.

    resetEnlargementFactor(enlargementFactor);
}














// Ellipsoid::~Ellipsoid()
//
// PURPOSE:
//      Class destructor.
//

Ellipsoid::~Ellipsoid()
{

}













// Ellipsoid::resetEnlargementFactor()
//
// PURPOSE: 
//      Compute covariance matrix, center coordinates, 
//      eigenvalues and eigenvectors matrix of the ellipsoid.
//      Eigenvalues are automatically enlarged according to the
//      enlargementFactor term.
//
// INPUT:
//      enlargementFactor: a double to contain the enlargement
//                         factor to be used for the chosen ellipsoid.
//
// OUTPUT:
//      void

void Ellipsoid::resetEnlargementFactor(const double newEnlargementFactor)
{
    // Save the new user-specified enlargement factor

    this->enlargementFactor = newEnlargementFactor;
    
    // Enlarge the eigenvalues with the user-specified factor

    enlargedEigenvalues = originalEigenvalues.sqrt() + enlargementFactor * originalEigenvalues.sqrt();
    enlargedEigenvalues = enlargedEigenvalues * enlargedEigenvalues;

    // Recompute the hypervolume contained in the enlarged ellipsoid

    hyperVolume = enlargedEigenvalues.prod();
    
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
//      ellipsoid: Ellipsoid object
//
// OUTPUT:
//      A boolean value specifying whether the two ellipsoids overlap (true) or not (false)
//
// REMARKS:
//      The coordinate system used by the ellipsoid in the argument must be the same as the
//      one used by this ellipsoid. This is relevant to know which coordinate corresponds to
//      which column in the covariance matrices of the ellipsoids.
//

bool Ellipsoid::overlapsWith(Ellipsoid ellipsoid)
{
    // Construct translation matrix

    MatrixXd T1 = MatrixXd::Identity(Ndimensions+1,Ndimensions+1);
    MatrixXd T2 = MatrixXd::Identity(Ndimensions+1,Ndimensions+1);
    
    T1.bottomLeftCorner(1,Ndimensions) = (-1.0) * this->getCenterCoordinates().transpose();
    T2.bottomLeftCorner(1,Ndimensions) = (-1.0) * ellipsoid.getCenterCoordinates().transpose();


    // Construct ellipsoid matrix in homogeneous coordinates

    MatrixXd A = MatrixXd::Zero(Ndimensions+1,Ndimensions+1);
    MatrixXd B = A;

    A(Ndimensions,Ndimensions) = -1;
    B(Ndimensions,Ndimensions) = -1;

    A.topLeftCorner(Ndimensions,Ndimensions) = this->getCovarianceMatrix().matrix().inverse();
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

    if (eigenSolver.info() != Success) abort();
    
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















// Ellipsoid::getSampleOfParameters()
//
// PURPOSE: 
//      Gets the protected data member sampleOfParameters.      
//
// OUTPUT:
//      An Eigen Array matrix of dimensions (Ndimensions, Nobjects) 
//      containing all the coordinates of the objects contained in the
//      ellipsoid.
//

ArrayXXd Ellipsoid::getSampleOfParameters()
{
    return sampleOfParameters;
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












// Ellipsoid::getNobjects()
//
// PURPOSE: 
//      Gets the protected data member Nobjects.      
//
// OUTPUT:
//      An integer containing the number of points (objects)
//      inside the ellipsoid.
//

int Ellipsoid::getNobjects()
{
    return Nobjects;
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













// Ellipsoid::getEnlargementFactor()
//
// PURPOSE: 
//      Get the private data member enlargementFactor.      
//
// OUTPUT:
//      A double containing the enlargement factor computed for the ellipsoid.
//

double Ellipsoid::getEnlargementFactor()
{
    return enlargementFactor;
}
