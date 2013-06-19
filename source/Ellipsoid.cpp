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
