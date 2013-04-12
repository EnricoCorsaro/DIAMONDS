#include "Ellipsoid.h"


// Ellipsoid::Ellipsoid()
//
// PURPOSE:
//      Class constructor.
//
// INPUT:
//      sampleOfParameters: an Eigen Array matrix of
//      size (Ndimensions, Nobjects), containing the
//      coordinates of the objects inside the ellipsoid.
//

Ellipsoid::Ellipsoid(RefArrayXXd sampleOfParameters, const int index)
: sampleOfParameters(sampleOfParameters),
  Nobjects(sampleOfParameters.cols()),
  index(index),
  Ndimensions(sampleOfParameters.rows())
{
    eigenvalues.resize(Ndimensions);
    centerCoordinates.resize(Ndimensions);
    eigenvectorsMatrix.resize(Ndimensions,Ndimensions);
    covarianceMatrix.resize(Ndimensions,Ndimensions);
}














// Ellipsoid::~Ellipsoid()
//
// PURPOSE:
//      Class destructor.
//

Ellipsoid::~Ellipsoid()
{

}













// Ellipsoid::getCenterCoordinates()
//
// PURPOSE: 
//      Compute covariance matrix, center coordinates, 
//      eigenvalues and eigenvectors matrix of the ellipsoid.
//      Eigenvalues are automatically enlarged according to the
//      enlargementFactor.
//
// OUTPUT:
//      void

void Ellipsoid::build(const double enlargementFactor)
{
    Functions::clusterCovariance(sampleOfParameters, covarianceMatrix, centerCoordinates);
    Functions::selfAdjointMatrixDecomposition(covarianceMatrix, eigenvalues, eigenvectorsMatrix);

    ArrayXd enlargedEigenvalues(Ndimensions);
    enlargedEigenvalues = eigenvalues.sqrt() + enlargementFactor*eigenvalues.sqrt();
    eigenvalues = enlargedEigenvalues * enlargedEigenvalues;
    hyperVolume = enlargedEigenvalues.prod();
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
//      An Eigen Array of dimensions (Ndimensions),
//      containing all original eigenvalues of the ellipsoid.
//

ArrayXd Ellipsoid::getEigenvalues()
{
    return eigenvalues;
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
//      Gets the protected data member covarianceMatrix.      
//
// OUTPUT:
//      An Eigen Array matrix of dimensions (Ndimensions, Ndimensions) 
//      containing the covariance matrix of the ellipsoid.
//

ArrayXXd Ellipsoid::getCovarianceMatrix()
{
    return covarianceMatrix;
}













// Ellipsoid::getEigenvectorsMatrix()
//
// PURPOSE: 
//      Gets the protected data member enlargedEigenvectors.      
//
// OUTPUT:
//      An Eigen Array matrix of dimensions (Ndimensions, Ndimensions), 
//      containing all eigenvectors of the ellipsoid.
//

ArrayXXd Ellipsoid::getEigenvectorsMatrix()
{
    return eigenvectorsMatrix;
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















// Ellipsoid::getIndex()
//
// PURPOSE: 
//      Gets the protected data member index.      
//
// OUTPUT:
//      An integer containing the reference number of the ellipsoid.
//

int Ellipsoid::getIndex()
{
    return index;
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
