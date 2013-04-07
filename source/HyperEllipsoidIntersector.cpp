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










// HyperEllipsoidIntersector::intersection()
//
// PURPOSE:
//      Determines whether two input ellipsoids are overlapping according to
//      the algorithm described by Alfano & Greer (2003; Journal of Guidance,
//      Control and Dynamics, 26, 1).
//
// INPUT:
//      covarianceMatrix1: ArrayXXd containing covariance matrix of the first ellipsoid
//      centerCoordinates1: ArrayXd containing center coordinates of first ellipsoid
//      covarianceMatrix2: ArrayXXd containing covariance matrix of the second ellipsoid
//      centerCoordinates2: ArrayXd containing center coordinates of second ellipsoid
//
// OUTPUT:
//      A boolean value specifying whether the two ellipsoids intersect (true) or not (false)
//
// REMARKS:
//      Coordinates of centers have to coincide in same order of the covariance matrix dimensions.
//      E.g. if first center coordinate is x, then first row and column in covariance matrix refer to x coordinate.
//

bool HyperEllipsoidIntersector::intersection(const RefArrayXXd covarianceMatrix1, const RefArrayXd centerCoordinates1, const RefArrayXXd covarianceMatrix2, const RefArrayXd centerCoordinates2)
{
    assert(covarianceMatrix1.cols() == covarianceMatrix2.cols());
    assert(centerCoordinates1.size() == centerCoordinates2.size());
    assert(covarianceMatrix1.cols() == centerCoordinates1.size());

    int Ndimensions = centerCoordinates1.size();


    // Construct translation matrix

    MatrixXd T1 = MatrixXd::Identity(Ndimensions+1,Ndimensions+1);
    MatrixXd T2 = MatrixXd::Identity(Ndimensions+1,Ndimensions+1);
    
    T1.bottomLeftCorner(1,Ndimensions) = -1.*centerCoordinates1.transpose();
    T2.bottomLeftCorner(1,Ndimensions) = -1.*centerCoordinates2.transpose();


    // Construct ellipsoid matrix in homogeneous coordinates

    MatrixXd A = MatrixXd::Zero(Ndimensions+1,Ndimensions+1);
    MatrixXd B = A;

    A(Ndimensions,Ndimensions) = -1;
    B(Ndimensions,Ndimensions) = -1;
    
    A.topLeftCorner(Ndimensions,Ndimensions) = covarianceMatrix1.matrix().inverse();
    B.topLeftCorner(Ndimensions,Ndimensions) = covarianceMatrix2.matrix().inverse();

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












// HyperEllipsoidIntersector::findNonOverlappingEllipsoids()
//
// PURPOSE:
//      Determines which of Nclusters input ellipsoids are not overlapping.
//
// INPUT:
//      Nclusters: an integer specifying the number of ellipsoids to check for
//      overlaps.
//      allEnlargedEigenvalues: an Eigen Array of dimensions (Nclusters * Ndimensions),
//      containing all enlarged eigenvalues of each ellipsoid.
//      allEigenvectorsMatrix: an Eigen Array matrix of dimensions 
//      (Ndimensions, Nclusters * Ndimensions), containing 
//      all eigenvectors of each ellipsoid.
//      allCentersCoordinates: An Eigen Array of dimensions (Nclusters*Ndimensions)
//      containing all the arrays of the center coordinates of each cluster.
//
// OUTPUT:
//      An Eigen Array of integers specifying the subscripts of the ellipsoids that 
//      are not overlapping. A single element array containing the value -1 is returned 
//      in case no non-overlapping ellipsoids are found.
//
ArrayXi HyperEllipsoidIntersector::findNonOverlappingEllipsoids(const int Nclusters, const RefArrayXd allEnlargedEigenvalues, const RefArrayXXd allEigenvectorsMatrix, const RefArrayXd allCentersCoordinates)
{
    assert(allEnlargedEigenvalues.size() == allCentersCoordinates.size());
    assert(allEigenvectorsMatrix.rows() == allCentersCoordinates.size()/Nclusters);

    int Ndimensions = allEigenvectorsMatrix.rows();

    ArrayXd enlargedEigenvalues1(Ndimensions);
    ArrayXd centerCoordinates1(Ndimensions);
    ArrayXXd eigenvectorsMatrix1(Ndimensions,Ndimensions);
    ArrayXXd covarianceMatrix1(Ndimensions,Ndimensions);
    ArrayXd enlargedEigenvalues2(Ndimensions);
    ArrayXd centerCoordinates2(Ndimensions);
    ArrayXXd eigenvectorsMatrix2(Ndimensions,Ndimensions);
    ArrayXXd covarianceMatrix2(Ndimensions,Ndimensions);
    
    bool overlap = false;
    bool noOverlapFlag;
    int count = 0;
    ArrayXi nonOverlappingEllipsoids;

    for (int i = 0; i < Nclusters-1; i++)
    {
        centerCoordinates1 = allCentersCoordinates.segment(i*Ndimensions, Ndimensions);
        enlargedEigenvalues1 = allEnlargedEigenvalues.segment(i*Ndimensions, Ndimensions);
        eigenvectorsMatrix1 = allEigenvectorsMatrix.block(0, i*Ndimensions, Ndimensions, Ndimensions);
        covarianceMatrix1 = eigenvectorsMatrix1.matrix() * enlargedEigenvalues1.matrix().asDiagonal() * eigenvectorsMatrix1.matrix().transpose();
        
        for (int j = i + 1; j < Nclusters; j++)
        {
            centerCoordinates2 = allCentersCoordinates.segment(j*Ndimensions, Ndimensions);
            enlargedEigenvalues2 = allEnlargedEigenvalues.segment(j*Ndimensions, Ndimensions);
            eigenvectorsMatrix2 = allEigenvectorsMatrix.block(0, j*Ndimensions, Ndimensions, Ndimensions);
            covarianceMatrix2 = eigenvectorsMatrix2.matrix() * enlargedEigenvalues2.matrix().asDiagonal() * eigenvectorsMatrix2.matrix().transpose();
            overlap = intersection(covarianceMatrix1, centerCoordinates1, covarianceMatrix2, centerCoordinates2);

            if (overlap)    // If overlap occurred, go to next i-th ellipsoid
            {
                noOverlapFlag = false;  // Set noOverlap flag to false
                break;
            }
            else    // If no overlap occurred, go to next j-th ellipsoid
            {
                noOverlapFlag = true;  // Set noOverlap flag to true
                continue;        
            }
                
        }
        
        if (noOverlapFlag)      // If no overlaps for the i-th ellipsoid are found...
        {
            if (i == Nclusters-2)       // save i and i+1 if i is the last one
            {
                nonOverlappingEllipsoids.conservativeResize(count+2);
                nonOverlappingEllipsoids(count) = i;
                nonOverlappingEllipsoids(count+1) = i+1;
                continue;
            }
            else                        // save i and go to next i otherwise
            {
                nonOverlappingEllipsoids.conservativeResize(count+1);
                nonOverlappingEllipsoids(count) = i;
                count++;
            }
        }
        else                    // If overlaps are found, go to next i
            continue;
    }

    return nonOverlappingEllipsoids;
}



















// HyperEllipsoidIntersector::checkPointForOverlap()
//
// PURPOSE:
//      Determines if the input point belongs to the input enlarged ellipsoid.
//
// INPUT:
//      enlargedEigenvalues: the enlarged eigenvalues of the ellipsoid.
//      eigenVectorsMatrix: the matrix of dimensions (Ndimensions, Ndimensions)
//      containing the column eigenvectors of the ellipsoid covariance matrix.
//      centersCoordinates: Eigen Array containing the center coordinates
//      of the input ellipsoid.
//      pointCoordinates: coordinates of the point to verify.
//
// OUTPUT:
//      A boolean value specifying whether the point belongs to the input
//      enlarged ellipsoid.
//

bool HyperEllipsoidIntersector::checkPointForOverlap(const RefArrayXd enlargedEigenvalues, const RefArrayXXd eigenVectorsMatrix, const RefArrayXd centerCoordinates, const RefArrayXd pointCoordinates)
{
    int Ndimensions = pointCoordinates.size();
       
    // Construct translation matrix

    MatrixXd T = MatrixXd::Identity(Ndimensions+1,Ndimensions+1);
    
    T.bottomLeftCorner(1,Ndimensions) = -1.*centerCoordinates.transpose();


    // Construct ellipsoid matrix in homogeneous coordinates

    MatrixXd A = MatrixXd::Zero(Ndimensions+1,Ndimensions+1);
    A(Ndimensions,Ndimensions) = -1;
    
    MatrixXd V = eigenVectorsMatrix.matrix();
    MatrixXd C = MatrixXd::Zero(Ndimensions, Ndimensions);
    
    C = V * enlargedEigenvalues.matrix().asDiagonal() * V.transpose();      // Covariance matrix
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
