#include "Functions.h"


// Functions::lorentzProfile()
//
// PURPOSE: 
//      Computes a simple Lorentzian profile given the centroid, the amplitude and the width.
//
// INPUT:
//      predictions : vector containing the result
//      covariates : vector containing independent variable values
//      centroid : centroid of the Lorentzian profile (default = 0)
//      amplitude : maximum height of the Lorentzian profile (default = 1)
//      gamma : width of the Lorentzian profile (default = 1)
//
// OUTPUT:
//      void
//
// REMARKS:
//      Saves the predictions into the input vector predictions.
//

void Functions::lorentzProfile(RefArrayXd predictions, const RefArrayXd covariates, 
                               const double centroid, const double amplitude, const double gamma)
{
    predictions = (amplitude*amplitude)/((covariates-centroid)*(covariates-centroid) + (gamma/2.)*(gamma/2.));
}











// Functions::modeProfile()
//
// PURPOSE: 
//      Computes a single Lorentzian profile that models an oscillation mode in a power spectrum.
//      The profile is computed given the centroid, the height and the linewidth.
//
// INPUT:
//      predictions : vector containing the result
//      covariates : vector containing independent variable values
//      centroid : centroid of the Lorentzian profile (default = 0)
//      height : maximum height of the Lorentzian profile (default = 1)
//      linewidth : width of the Lorentzian profile (mode linewidth) (default = 1)
//
// OUTPUT:
//      void
//
// REMARKS:
//      Saves the predictions into the input vector predictions.
//

void Functions::modeProfile(RefArrayXd predictions, const RefArrayXd covariates, 
                               const double centroid, const double height, const double linewidth)
{
    predictions = height/(1.0 + (4.0*(covariates-centroid).square()/(linewidth*linewidth)));
}











// Functions::logGaussProfile()
//
// PURPOSE: 
//      Computes the logarithm of a Gaussian profile given the centroid, 
//      the standard deviation and the amplitude, for one given x-value
//
// INPUT:
//      covariate : independent variable
//      mu : mean value of the Gaussian profile (default = 0)
//      sigma : standard deviation of the Gaussian profile (default = 1)
//      amplitude : maximum amplitude of the Gaussian profile (default = 1)
//
// OUTPUT:
//      The natural logarithm of the Gaussian profile value of the independent variable.
//

double Functions::logGaussProfile(const double covariate, const double mu, 
                                  const double sigma, const double amplitude)
{
    const double prefactor = log(amplitude) - 0.5 * log(2*Functions::PI) - log(sigma);

    return prefactor - 0.5 * (covariate - mu) * (covariate - mu) / (sigma * sigma);

} 








// Functions::logGaussProfile()
//
// PURPOSE: 
//      Computes the logarithm of a Gaussian profile given the centroid, 
//      the standard deviation and the amplitude for a set of values (overloaded)
//
// INPUT:
//      predictions : Eigen array containing the result
//      covariates : Eigen array containing independent variable values
//      mu : mean value of the Gaussian profile (default = 0)
//      sigma : standard deviation of the Gaussian profile (default = 1)
//      amplitude : maximum amplitude of the Gaussian profile (default = 1)
//
// OUTPUT:
//      void
//

void Functions::logGaussProfile(RefArrayXd predictions, const RefArrayXd covariates, const double mu, 
                                const double sigma, const double amplitude)
{
    const double prefactor = log(amplitude) - 0.5 * log(2*Functions::PI) - log(sigma);
    predictions = prefactor - 0.5 * (covariates - mu) * (covariates - mu) / (sigma * sigma);
    
}









// Functions::logGaussLikelihood() 
//
// PURPOSE: 
//      Computes the log-Gaussian likelihood from a set of observations,
//      uncertainties, and theoretical predictions.
//
// INPUT:
//      observations : array containing the oberved values
//      predictions : array containing the corresponding predicted values of x
//      uncertainties : array containing the uncertanties on the observed values
//
// OUTPUT:
//      The resulting value of the log-Gaussian likelihood.
//

double Functions::logGaussLikelihood(const RefArrayXd observations, const RefArrayXd predictions, const RefArrayXd uncertainties)
{
    if ((observations.size() != predictions.size()) || (observations.size() != uncertainties.size()))
    {
        cerr << "Array dimensions do not match. Quitting program." << endl;
        exit(EXIT_FAILURE);
    }
    
    ArrayXd delta;
    ArrayXd lambda0;
    ArrayXd lambda;
    
    delta = ((observations - predictions)*(observations - predictions)) / (uncertainties*uncertainties);
    lambda0 = -1.*log(sqrt(2.*PI) * uncertainties);
    lambda = lambda0 -0.5*delta;
    
    return lambda.sum();

}










// Functions::clusterCovariance()
//
// PURPOSE:
//      Compute the covariance matrix of a set of points belonging to a cluster.
//
// INPUT:
//      clusterSample: an Eigen Array matrix of size (Ndimensions, Npoints)
//      containing coordinates of all the points belonging to the cluster.
//      covarianceMatrix: an Eigen Array matrix of size (Ndimensions, Ndimensions) where the
//      covariance matrix is stored.
//      centerCoordinates: an Eigen Array of size Ndimensions where the coordinates of the
//      center (mean values) of the cluster of points are stored.
//
// OUTPUT:
//      void
//

void Functions::clusterCovariance(RefArrayXXd const clusterSample, RefArrayXXd covarianceMatrix, 
                                  RefArrayXd centerCoordinates)
{    
    int Ndimensions = clusterSample.rows();
    int Npoints = clusterSample.cols();
    covarianceMatrix.resize(Ndimensions, Ndimensions);
    centerCoordinates.resize(Ndimensions);

    for (int i=0; i < Ndimensions; i++)
    {
        centerCoordinates(i) = (clusterSample.row(i).sum())/Npoints;
    }

    double biasFactor = 1./(Npoints-1);

    for (int i=0; i < Ndimensions; i++)
    {
        for (int j=0; j < Ndimensions; j++)
        {
            covarianceMatrix(i,j) = biasFactor * ((clusterSample.row(i) - centerCoordinates(i))*(clusterSample.row(j) - centerCoordinates(j))).sum();
        }
    }

}











// Functions::selfAdjointMatrixDecomposition()
//
// PURPOSE:
//      Compute the decomposition of a covariance matrix into eigenvectors and eigenvalues.
//
// INPUT:
//      covarianceMatrix:       an Eigen Array matrix to be decomposed
//                              eigenvalues: an Eigen Array to contain the eigenvalues 
//                              of the covariance matrix.
//      eigenvalues:            an Eigen Array containing the eigenvalues of the covariance matrix
//                              listed in ascending order.
//      eigenvectorsMatrix:     an Eigen Array matrix to contain the matrix of
//                              the eigenvectors of the covariance matrix.
//
// OUTPUT:
//      A boolean value that is true if the decomposition was successfull, false otherwise.
//

bool Functions::selfAdjointMatrixDecomposition(RefArrayXXd const covarianceMatrix, RefArrayXd eigenvalues, 
                                               RefArrayXXd eigenvectorsMatrix)
{
    assert(covarianceMatrix.cols() == covarianceMatrix.rows());
    assert(eigenvalues.size() == covarianceMatrix.cols());
    assert(eigenvectorsMatrix.cols() == eigenvectorsMatrix.rows());
    assert(eigenvectorsMatrix.cols() == eigenvalues.size());

    SelfAdjointEigenSolver<MatrixXd> eigenSolver(covarianceMatrix.matrix());

    if (eigenSolver.info() != Success)
    {
        cout << "Covariance Matrix decomposition failed." << endl;
        cout << "Quitting program" << endl;
        return false;
    }

    eigenvalues = eigenSolver.eigenvalues();
    eigenvectorsMatrix = eigenSolver.eigenvectors();

    return true;
}










// Functions::product()
//
// PURPOSE: 
//      Computes the product of the elements contained in a vector of doubles of class vector.
//
// INPUT:
//      vec : vector of values to be multiplied
//
// OUTPUT:
//      The productoria of the vector elements.
//

inline double Functions::product(const vector<double> &vec)
{
    return accumulate(vec.begin(), vec.end(), 1.0, multiplies<double>());
}









// Functions::sum()
//
// PURPOSE: 
//      Computes the sum of the elements contained in a vector of doubles of vector class.
//
// INPUT:
//      vec : vector of values to be added.
//
// OUTPUT:
//      The summation of the vector elements.
//

inline double Functions::sum(const vector<double> &vec)
{
    return accumulate(vec.begin(), vec.end(), 0.0);
}










// Functions::logExpSum()
//
// PURPOSE: 
//      Computes a logaritmic summation of exponentials in order to avoid overflow errors.
//
// INPUT: 
//      x : a first variable to be added
//      y : a second variable to be added
//
// OUTPUT: 
//      The logarithmic summation of the exponentials log(exp(x)+exp(y)).
//

double Functions::logExpSum(const double x, const double y)
{
    return (x >= y ? x + log(1.+exp(y-x)) : y + log(1.+exp(x-y)));
}









// Functions::logExpDifference()
//
// PURPOSE: 
//      Computes a logaritmic difference of exponentials in order to avoid overflow errors.
//
// INPUT: 
//      x : a first variable to be added
//      y : a second variable to be added
//
// OUTPUT: 
//      The logarithmic difference of the exponentials log(exp(x)-exp(y)).
//

double Functions::logExpDifference(const double x, const double y)
{
    return (x >= y ? x + log(1.0 - exp(y-x)) : y + log(1.0 - exp(x-y)));
}












// Functions::sortElementsDouble()
//
// PURPOSE: 
//      Sorts the element of the first input array in increasing order 
//      and the elements of the second input array according to the 
//      sorting of the first array.
//
// INPUT: 
//      array1: a first Eigen Array of double numbers to be sorted in increasing order 
//      array2: a second Eigen Array to be sorted according to array1
//
// OUTPUT: 
//      void
//

void Functions::sortElementsDouble(RefArrayXd array1, RefArrayXd array2)
{
    for (int i = 0; i < array1.size(); i++)
    {
        for (int j = 1; j < (array1.size()-i); j++)
        {
            if (array1(j-1) > array1(j))
            {
                SWAPDOUBLE(array1(j-1),array1(j));        // SWAP array1 elements in increasing order
                SWAPDOUBLE(array2(j-1),array2(j));        // SWAP array2 elements accordingly
            }
            else
                if (array1(j-1) == array1(j))
                    continue;
        }
    }
    
}













// Functions::sortElementsInt()
//
// PURPOSE: 
//      Sorts the element of the first input array in increasing order 
//      and the elements of the second input array according to the 
//      sorting of the first array.
//
// INPUT: 
//      array1: vector of integer values to be sorted in increasing order 
//      array2: a second Eigen Array to be sorted according to array1
//
// OUTPUT: 
//      void
//

void Functions::sortElementsInt(vector<int> &array1, RefArrayXd array2)
{
    for (int i = 0; i < array1.size(); i++)
    {
        for (int j = 1; j < (array1.size()-i); j++)
        {
            if (array1[j-1] > array1[j])
            {
                SWAPINT(array1[j-1], array1[j]);         // SWAP array1 elements in increasing order
                SWAPDOUBLE(array2(j-1), array2(j));      // SWAP array2 elements accordingly
            }
            else
                if (array1[j-1] == array1[j])
                    continue;
        }
    }
    
}












// Functions::findArrayIndicesWithinBoundaries()
//
// PURPOSE: 
//      This function returns an array containing the indices of the input array that correspond
//      to the elemtns contained within the input boundaries.
//
// INPUT:
//      array:       an Eigen array whose indices we want to find
//      lowerBound:  a double specifying the smallest value allowed in the search
//      upperBound:  a double specifying the largest value allowed in the search
//
// OUTPUT:
//      a vector<int> containing the indices of the input array that correspond to its elements
//      that fall within the input boundaries.
//
// REMARKS:
//      There is no requirement for the input array to have its element sorted in any order. 
//      The finding of the indices is done correctly independently of the sorting of the input array elements.
//

vector<int> Functions::findArrayIndicesWithinBoundaries(RefArrayXd const array, double lowerBound, double upperBound)
{
    // At least one point is needed

    assert(array.size() >= 1);
    vector<int> arrayIndicesWithinBoundaries;

    for (int i = 0; i < array.size(); i++)
    {
        if ((array(i) >= lowerBound) && (array(i) <= upperBound))
        {
            arrayIndicesWithinBoundaries.push_back(i);
        }
        else
            continue;
    }

    return arrayIndicesWithinBoundaries;
}










// Functions::AkimaInterpolation()
// 
// PURPOSE:
//      This function computes an Akima cubic spline interpolation of some one-dimensional input data,
//      given a grid of new values we want the input data to be interpolated. The algorithm is documented
//      by Martin Rottinger, The Akima Interpolation, 1999 and by David Eberly, Akima Interpolation 
//      for Nonuniform 1D Data, 2013.
//      The advantage of Akima interpolation is that it is less affected by outliers in the input data
//      and it is able to reproduce a more natural interpolation, having less wiggling with respect to the standard
//      cubic spline interpolation.
//
// INPUT:
//      observedAbscissa:       an Eigen array containing the input data covariates to be used for the interpolation
//      observedOrdinate:       an Eigen array containing the corresponding dependent variables of the input abscissa
//      interpolatedAbscissa:   an Eigen array containing the values of the covariates for which we need to compute the
//                              interpolated dependent variables
//
// OUTPUT:
//      An eigen array of doubles containing the newly compute ordinates (dependent variables) of the corresponding
//      input grid of covariates for which the interpolation was required.
//
// REMARKS:
//      It is required that all input abscissa arrays contain values sorted in increasing order.
//      Obivously, input ordinates have to be sorted accordingly to their corresponding abscissa.
//      In addition, observed abscissa boundaries must comprise those of the interpolated abscissa.
//      This means that no extrapolation can be done.
//

ArrayXd Functions::AkimaInterpolation(RefArrayXd const observedAbscissa, RefArrayXd const observedOrdinate, RefArrayXd const interpolatedAbscissa)
{
    // Number of data points. The number of intervals is clearly size-1
    
    int size = observedAbscissa.size();           
    
    
    // Number of interpolation grid points.
    
    int interpolatedSize = interpolatedAbscissa.size();
 

    // Since the formula requires at least 5 data points, check if array size is not below 5,
    // if the interpolated grid has at least one point, and if input abscissa and ordinate 
    // have same number of elements.
    
    assert(size >= 5);
    assert(interpolatedSize >= 1);

    
    // Check if boundaries set by observed data points are respected by interpolated points
    
    assert (observedAbscissa(0) <= interpolatedAbscissa(0));
    assert (observedAbscissa(size-1) >= interpolatedAbscissa(size-1)); 


    // Initialize arrays of differences in both ordinate and abscissa

    ArrayXd differenceOrdinate = observedOrdinate.segment(1,size-1) - observedOrdinate.segment(0,size-1);
    ArrayXd differenceAbscissa = observedAbscissa.segment(1,size-1) - observedAbscissa.segment(0,size-1);
    

    // Compute the ratios for all the input values. These ratios are in number "size-1" (from ratio_2 to ratio_(size))
    // + 2 ratios to the left (ratio_0, ratio_1) + 2 ratio to the right (ratio_(size+1), ratio_(size+2)). 
    // The total number of ratios is then size+3. 
    
    int ratiosSize = size + 3;
    ArrayXd ratios = ArrayXd::Zero(ratiosSize);
    ratios.segment(2,size-1) = differenceOrdinate/differenceAbscissa;
    ratios(1) = (2 * ratios(2)) - ratios(3);                                         // ratio_1
    ratios(0) = (2 * ratios(1)) - ratios(2);                                         // ratio_0
    ratios(ratiosSize-2) = (2 * ratios(ratiosSize-3)) - ratios(ratiosSize-4);        // ratio_(size+1)
    ratios(ratiosSize-1) = (2 * ratios(ratiosSize-2)) - ratios(ratiosSize-3);        // ratio_(size+2)

    
    // Compute weights to be used in the formula of the first derivatives. 
    // For each derivative, two weights are required.

    ArrayXd weightsLeft(size);
    ArrayXd weightsRight(size);
    weightsLeft = (ratios.segment(3,size) - ratios.segment(2,size)).abs();
    weightsRight = (ratios.segment(1,size) - ratios.segment(0,size)).abs();


    // Compute the first derivatives at each data point

    vector<double> firstDerivatives(size);

    for (int i = 0; i < size; i++)
    {
        if ((weightsLeft(i) == weightsRight(i)) == 0)
        {
            // Both weights are zero, hence adopt the arithmetic average

            firstDerivatives[i] = 0.5 * ( ratios(i+1) + ratios(i+2) );
        }
        else
        {
            // Since at least one of the weights is != 0, adopt the Akima rule for the first derivative

            firstDerivatives[i] = (weightsLeft(i)*ratios(i+1) + weightsRight(i)*ratios(i+2)) /
                                  (weightsLeft(i) + weightsRight(i));
        }
    }
    

    // Start a loop over each bin of the input data point, i.e. between 
    // observedOrdinate(i) and observedOrdinate(i+1)

    ArrayXd interpolatedOrdinate(interpolatedSize);
    ArrayXd remainingInterpolatedAbscissa = interpolatedAbscissa;       // The remaining part of the array of interpolated abscissa
    int cumulatedBinSize = 0;               // The cumulated number of interpolated points from the beginning

    for (int i = 0; i < size-1; i++)
    {
        // For each bin given by the input data points array (we have size-1 bins in total), 
        // compute the coefficients of the corresponding cubic polynomial and 
        // the new ordinates for the corresponding interpolated abscissa.
        
        double coeff0 = observedOrdinate(i);
        double coeff1 = firstDerivatives[i];
        double coeff2 = 3*(observedOrdinate(i+1) - observedOrdinate(i) - firstDerivatives[i] ) - 
                        (firstDerivatives[i+1] - firstDerivatives[i]);
        double coeff3 = (firstDerivatives[i+1] - firstDerivatives[i]) - 2*(observedOrdinate(i+1) - 
                        observedOrdinate(i) - firstDerivatives[i]);


        // Find which values of interpolatedAbscissa are containined within the selected bin of observedAbscissa.
        // Since elements in interpolatedAbscissa are monotonically increasing, we cut the input array each time
        // we identify the elements of the current bin. This allows to speed up the process.

        double lowerBound = observedAbscissa(i);
        double upperBound = observedAbscissa(i+1);
        vector<int> selectedIndices = findArrayIndicesWithinBoundaries(remainingInterpolatedAbscissa, lowerBound, upperBound);

        int binSize = selectedIndices.size();       // Total number of interpolated points falling in the current bin
        double binLength = differenceAbscissa(i);
        ArrayXd interpolatedAbscissaInCurrentBin = remainingInterpolatedAbscissa.segment(0, binSize);

        
        // Compute interpolated ordinates for the current bin

        ArrayXd normalizedAbscissa = (interpolatedAbscissaInCurrentBin - lowerBound)/binLength;
        ArrayXd interpolatedOrdinateInCurrentBin = coeff0 + coeff1*normalizedAbscissa + coeff2*normalizedAbscissa.square() +
                                                   coeff3*normalizedAbscissa.cube();
       
        
        // Merge bin ordinate into total array of ordinate

        interpolatedOrdinate.segment(cumulatedBinSize, binSize) = interpolatedOrdinateInCurrentBin;
        

        // Reduce size of array remainingInterpolatedAbscissa by binSize elements and initialize the array
        // with remaining part of interpolatedAbscissa
        
        int currentRemainingSize = remainingInterpolatedAbscissa.size();
        remainingInterpolatedAbscissa.resize(currentRemainingSize - binSize);
        cumulatedBinSize += binSize;
        remainingInterpolatedAbscissa = interpolatedAbscissa.segment(cumulatedBinSize, interpolatedSize-cumulatedBinSize);
    }

    return interpolatedOrdinate;
}
