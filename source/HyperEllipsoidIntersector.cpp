#include "HyperEllipsoidIntersector.h"

// HyperEllipsoidIntersector::HyperEllipsoidIntersector()
//
// PURPOSE: 
//      Derived class constructor.
//

HyperEllipsoidIntersector::HyperEllipsoidIntersector()
{

}











// HyperEllipsoidIntersector::~HyperEllipsoidIntersector()
//
// PURPOSE: 
//      Derived class destructor.
//

HyperEllipsoidIntersector::~HyperEllipsoidIntersector()
{

}










// HyperEllipsoidIntersector::result()
//
// PURPOSE:
//      Determines whether two input ellipsoids are overlapping according to
//      the algorithm described by Alfano & Greer (2003; Journal of Guidance,
//      Control and Dynamics, 26, 1).
//
// INPUT:
//      covarianceMatrix1: ArrayXXd containing covariance matrix of the first ellipsoid
//      covarianceMatrix2: ArrayXXd containing covariance matrix of the second ellipsoid
//      centerCoordinates1: ArrayXd containing center coordinates of first ellipsoid
//      centerCoordinates2: ArrayXd containing center coordinates of second ellipsoid
//
// OUTPUT:
//      An integer specifying whether the two ellipsoids intersect (1) or not (0)
//
// REMARKS:
//      Coordinates of centers have to coincide in same order of the covariance matrix dimensions.
//      E.g. if first center coordinate is x, then first row and column in covariance matrix refer to x coordinate.
//

int HyperEllipsoidIntersector::result(RefArrayXXd covarianceMatrix1, RefArrayXXd covarianceMatrix2, RefArrayXd centerCoordinates1, RefArrayXd centerCoordinates2)
{
    assert(covarianceMatrix1.cols() == covarianceMatrix2.cols());
    assert(centerCoordinates1.size() == centerCoordinates2.size());
    assert(covarianceMatrix1.cols() == centerCoordinates1.size());

    int Ndim = centerCoordinates1.size()


    // Construct translation matrix

    MatrixXd T1 = MatrixXd::Identity(Ndim+1,Ndim+1);
    MatrixXd T2 = MatrixXd::Identity(Ndim+1,Ndim+1);
    
    T1.bottomLeftCorner(1,Ndim) = -1.*centerCoordinates1.transpose();
    T2.bottomLeftCorner(1,Ndim) = -1.*centerCoordinates2.transpose();


    // Construct ellipsoid matrix in homogeneous coordinates

    MatrixXd A = MatrixXd::Zero(Ndim+1,Ndim+1);
    MatrixXd B = A;

    A(Ndim,Ndim) = -1;
    B(Ndim,Ndim) = -1;
    
    A.topLeftCorner(Ndim,Ndim) = covarianceMatrix1.matrix().inverse();
    B.topLeftCorner(Ndim,Ndim) = covarianceMatrix2.matrix().inverse();

    MatrixXd AT = T*A*T.transpose();        // Translating to ellispoid center
    MatrixXd BT = T*B*T.transpose();        // Translating to ellispoid center


    // Compute Hyper Quadric Matrix generating from the two ellipsoids 
    // and derive its eigenvalues decomposition

    MatrixXd C = AT.inverse() * BT;
    MatrixXcd CC(Ndim+1,Ndim+1);

    CC.imag() = MatrixXd::Zero(Ndim+1,Ndim+1); 
    CC.real() = C;
    
    ComplexEigenSolver<MatrixXcd> eigenSolver(CC);

    if (eigensolver.info() != Success) abort();
    
    MatrixXcd E = eigenSolver.eigenvalues();
    MatrixXcd V = eigenSolver.eigenvectors();

    int intersection = 0;       // Start with no intersection
    double pointA;              // Point laying in elliposid A
    double pointB;              // Point laying in ellipsoid B
    
    for (int i = 0; i < Ndim+1; i++)      // Loop over all eigenvectors
    {
        if (V(Ndim,i).real() == 0)      // Skip inadmissible eigenvectors
            continue;                   
        else if (E(i).imag() != 0)
            {
                V.col(i) = V.col(i).array() * (V.conjugate())(Ndim,i);      // Multiply eigenvector by complex conjugate of last element
                V.col(i) = V.col(i).array() / V(Ndim,i).real();             // Normalize eigenvector to last component value
                pointA = V.col(i).transpose().real() * AT * V.col(i).real();        // Evaluate point from elliposid A
                pointB = V.col(i).transpose().real() * BT * V.col(i).real();        // Evaluate point from ellipsoid B

                if ((pointA <= 0) && (pointB <= 0))     // Accept only if point belongs to both ellipsoids
                {
                    intersection = 1;                   // Exit if intersection is found
                    break;
                }
            }
    }

    return intersection;

}




