#include "PrincipalComponentProjector.h"

// PrincipalComponentProjector:PrincipalComponentProjector()
//
// PURPOSE: 
//      Class constructor for the derived class PrincipalComponent of the Projector abstract based class.
//
// INPUT:
//      printNdimensions: a boolean specifying whether one wants the reduced number
//      of dimensions to be explicitly printed on the screen.
//      scalingFactor: a double specifying the value of the scaling of the average eigenvalue (according to the correction
//      applied by Jolliffe to the Kaiser-Guttman test. A value of 0.7 is set by default.
//

PrincipalComponentProjector::PrincipalComponentProjector(bool printNdimensions, double scalingFactor)
: Projector(printNdimensions),
  scalingFactor(scalingFactor)
{

}











// PrincipalComponentProjector::projection()
//
// PURPOSE: 
//      Given a sample of N-dimensional points, use the linear Principal Component Analysis
//      to identify the principal components (PCs) of the original dataset. Subsequently select
//      only those component that contain the largest amount of variance in the data through the
//      Kaiser-Guttman test, modified by Jolliffe (1972) to consider PCs having eigenvalues > 70% of the average eigenvalue.
//      
//
// INPUT:
//      sample(Ndimensions, Npoints): sample of N-dimensional points
// 
// OUTPUT:
//      The reduced k-dimensional sample, with k <= N, but having the same number of points.
//

ArrayXXd PrincipalComponentProjector::projection(RefArrayXXd sample)
{
    Ndimensions = sample.rows();
    Npoints = sample.cols();
  
    
    // Evaluate the center of the sample (sample mean) in each dimensions
    
    ArrayXd sampleMean(Ndimensions);
    sampleMean = sample.rowwise().mean();       // Column array containing mean values for each dimension
    
    ArrayXXd differenceFromCenters(Ndimensions, Npoints);
    differenceFromCenters = sample.colwise() - sampleMean;

    
    // Evaluate the covariance matrix of the entire sample (N-dimensional matrix)
    
    double biasFactor = 1.0/(Npoints-1.0);
    ArrayXXd covarianceMatrix(Ndimensions, Ndimensions);
    covarianceMatrix = differenceFromCenters.matrix()*differenceFromCenters.matrix().transpose() * biasFactor;

    ArrayXXd varianceArray(Ndimensions,1);
    varianceArray = covarianceMatrix.matrix().diagonal();

    
    // Evaluate the correlation matrix of the entire sample so that its trace equals the number of dimensions, hence the variances
    // are normalized by the total variance from each dimension

    ArrayXXd correlationMatrix(Ndimensions, Ndimensions);
    correlationMatrix = varianceArray.sqrt().matrix().asDiagonal().inverse() * covarianceMatrix.matrix() * 
                        varianceArray.sqrt().matrix().asDiagonal().inverse();
    

    // Compute the eigenvalues and eigenvectors of the correlation matrix
    
    ArrayXXd eigenvectors(Ndimensions, Ndimensions);
    ArrayXd eigenvalues(Ndimensions);
    
    
    // Extracted eigenvectors are column vectors and eigenvalues are in increasing order
    
    bool correlationMatrixDecompositionIsSuccessful = Functions::selfAdjointMatrixDecomposition(correlationMatrix, eigenvalues, eigenvectors);

    if (!correlationMatrixDecompositionIsSuccessful)
    {   
        cerr << "Error in decomposition from PCA correlation matrix." << endl;
        abort();
    }
    

    // Sort eigenvalues in decreasing order, and eigenvectors accordingly

    ArrayXXd zetaScores(Ndimensions, Npoints);
    zetaScores = eigenvectors.matrix().transpose() * sample.matrix();


    // Compute predictions for the given dimension according to the Kaiser-Guttman test modified by Jolliffe

    double scaledAverageEigenvalue = scalingFactor*eigenvalues.mean();
    reducedNdimensions = ((eigenvalues - scaledAverageEigenvalue) > 0).count();
    ArrayXXd reducedSample(reducedNdimensions, Npoints);
  

    // Consider a subset of the original sample with a number of features (dimensions) equal to reducedNdimensions,
    // and take the zeta scores corresponding to the largest eigenvalues found, associated to the retained PCs

    reducedSample = zetaScores.block(Ndimensions-reducedNdimensions,0,reducedNdimensions,Npoints);

    if (printNdimensions)
    {
        cout << "Number of effective dimensions from PCA: " << reducedNdimensions << endl;
    }

    // Return the reduced sample with scores orderder by decreasing eigenvalue

    return reducedSample.colwise().reverse();
}
