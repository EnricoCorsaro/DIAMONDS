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
//      centroid : centroid of the Lorentzian profile (default = 0), expressed in muHz
//      height : height of the Lorentzian profile (default = 1), expressed in ppm^2 / muHz
//      linewidth : width of the Lorentzian profile (mode linewidth) (default = 1), expressed in muHz
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











// Functions::modeProfileWithAmplitude()
//
// PURPOSE: 
//      Computes a single Lorentzian profile that models an oscillation mode in a power spectrum.
//      The profile is computed given the centroid, the amplitude and the linewidth.
//
// INPUT:
//      predictions : vector containing the result
//      covariates : vector containing independent variable values
//      centroid : centroid of the Lorentzian profile (default = 0), expressed in muHz
//      amplitude : amplitude of the Lorentzian profile (default = 1), expressed in ppm
//      linewidth : width of the Lorentzian profile (mode linewidth) (default = 1), expressed in muHz
//
// OUTPUT:
//      void
//
// REMARKS:
//      Saves the predictions into the input vector predictions.
//

void Functions::modeProfileWithAmplitude(RefArrayXd predictions, const RefArrayXd covariates, 
                               const double centroid, const double amplitude, const double linewidth)
{
    predictions = amplitude*amplitude/(Functions::PI * linewidth)/(1.0 + (4.0*(covariates-centroid).square()/(linewidth*linewidth)));
}












// Functions::modeProfileSinc()
//
// PURPOSE: 
//      Computes a normalized single sinc-square profile that models an unresolved oscillation mode in the power spectrum.
//      The profile is computed given the centroid, the height and the resolution of the dataset.
//
// INPUT:
//      predictions : vector containing the result
//      covariates : vector containing independent variable values
//      centroid : centroid of the Lorentzian profile (default = 0), expressed in muHz
//      height : height of the sinc-square profile (default = 1), expressed in ppm^2 / muHz
//      resolution : the distance from the centroid to the first zero of the sinc-square (default = 1), expressed in muHz
//
// OUTPUT:
//      void
//
// REMARKS:
//      Saves the predictions into the input vector predictions.
//

void Functions::modeProfileSinc(RefArrayXd predictions, const RefArrayXd covariates, 
                               const double centroid, const double height, const double resolution)
{
    ArrayXd sincFunctionArgument = Functions::PI*(covariates - centroid)/resolution;
    ArrayXd sincFunction = sincFunctionArgument.sin() / sincFunctionArgument;


    // Multiply the profile by the height in the PSD

    predictions = height*sincFunction.square();
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














// Functions::topDownMerge()
//
// PURPOSE: 
//      This function does the sorting of the input array1 and sorts the elements of
//      array2 according to this sorting.
//
// INPUT: 
//      array1: a first Eigen Array of double numbers to be sorted in increasing order 
//      array2: a second Eigen Array to be sorted according to array1
//      arrayCopy1: a copy of array1 containing the effect of the merge sorting
//      arrayCopy2: a copy of array2 reflecting the sorting performed on arrayCopy1
//      beginIndex: an integer specifying the starting index of the sorting process
//      middleIndex: an integer specifying the mid-point of the array to be sorted
//      endIndex: an integer specifying the ending index of the sorting process
//
// OUTPUT: 
//      void
//
// REMARK:
//      This function is not intented to be used separately
//      but only through the call to topDownMergeSort. It is used
//      within the function topDownSplitMerge.
//

void Functions::topDownMerge(RefArrayXd array1, RefArrayXd arrayCopy1, RefArrayXd array2, RefArrayXd arrayCopy2, int beginIndex, int middleIndex, int endIndex)
{
    // Specify the initial indices of first and second-half of the input array

    int firstIndexFirstPart = beginIndex;
    int firstIndexSecondPart = middleIndex;


    // Do the sorting

    for (int i = beginIndex; i < endIndex; i++)
    {
        if ((firstIndexFirstPart < middleIndex) && (firstIndexSecondPart >= endIndex || (array1(firstIndexFirstPart) <= array1(firstIndexSecondPart))))
        {
            arrayCopy1(i) = array1(firstIndexFirstPart);
            arrayCopy2(i) = array2(firstIndexFirstPart);
            firstIndexFirstPart++;      // Move to the next element
        }
        else
        {
            arrayCopy1(i) = array1(firstIndexSecondPart);
            arrayCopy2(i) = array2(firstIndexSecondPart);
            firstIndexSecondPart++;     // Move to the next element
        }
    }
}












// Functions::topDownMergeSort()
//
// PURPOSE: 
//      This function splits the input array in two parts and repeat the process 
//      until only 1-element segments are found. Then sorts all the segments and
//      merges them into a total, sorted array.
//
// INPUT: 
//      array1: a first Eigen Array of double numbers to be sorted in increasing order 
//      array2: a second Eigen Array to be sorted according to array1
//      arrayCopy1: a copy of array1 containing the effect of the merge sorting
//      arrayCopy2: a copy of array2 reflecting the sorting performed on arrayCopy1
//      beginIndex: an integer specifying the starting index of the sorting process
//      endIndex: an integer specifying the ending index of the sorting process
//
// OUTPUT: 
//      void
//
// REMARK:
//      This function calls itself recursively. It is not intented to be used separately
//      but only through the call to topDownMergeSort.
//

void Functions::topDownSplitMerge(RefArrayXd array1, RefArrayXd arrayCopy1, RefArrayXd array2, RefArrayXd arrayCopy2, int beginIndex, int endIndex)
{
    // If input array1 contains only 1 element it is already sorted

    if (endIndex - beginIndex < 2)
        return;


    // Find mid-point of the input array

    int middleIndex = (endIndex + beginIndex) / 2;
    
    
    // Separate first-half and second-half of the array and repeat the process

    topDownSplitMerge(array1, arrayCopy1, array2, arrayCopy2, beginIndex, middleIndex);
    topDownSplitMerge(array1, arrayCopy1, array2, arrayCopy2, middleIndex, endIndex);

    
    // Order elements in each half and merge them into a single, sorted, array

    topDownMerge(array1, arrayCopy1, array2, arrayCopy2, beginIndex, middleIndex, endIndex);
   

    // Copy elements of array sorted into original array
    
    int length = endIndex - beginIndex;
    array1.segment(beginIndex, length) = arrayCopy1.segment(beginIndex, length);
    array2.segment(beginIndex, length) = arrayCopy2.segment(beginIndex, length);
}












// Functions::topDownMergeSort()
//
// PURPOSE: 
//      Sorts the element of the first input array in increasing order 
//      and the elements of the second input array according to the 
//      sorting of the first array. This is done according to the
//      top-down implementation of the mergesort algorithm.
//
// INPUT: 
//      array1: a first Eigen Array of double numbers to be sorted in increasing order 
//      array2: a second Eigen Array to be sorted according to array1
//
// OUTPUT: 
//      void
//
// REMARK:
//      This function uses the separate function topDownSplitMerge defined above.
//

void Functions::topDownMergeSort(RefArrayXd array1, RefArrayXd array2)
{
    assert(array1.size() == array2.size());
    ArrayXd arrayCopy1 = array1;
    ArrayXd arrayCopy2 = array2;
    topDownSplitMerge(array1, arrayCopy1, array2, arrayCopy2, 0, arrayCopy1.size());
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
    assert(array1.size() == array2.size());
    
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
//      to the elements contained within the input boundaries.
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

    if (lowerBound < upperBound)
    {
        for (int i = 0; i < array.size(); ++i)
        {
            if ((array(i) >= lowerBound) && (array(i) <= upperBound))
                arrayIndicesWithinBoundaries.push_back(i);
            else
                continue;
        }
    }

    return arrayIndicesWithinBoundaries;
}












// Functions::countArrayIndicesWithinBoundaries()
//
// PURPOSE: 
//      This function counts the number of elements of an input array that
//      fall within the input boundaries.
//
// INPUT:
//      array:       an Eigen array whose indices we want to find
//      lowerBound:  a double specifying the smallest value allowed in the search
//      upperBound:  a double specifying the largest value allowed in the search
//
// OUTPUT:
//      an integer containing the indices of the input array that correspond to its elements
//      that fall within the input boundaries.
//
// REMARKS:
//      There is no requirement for the input array to have its element sorted in any order. 
//      The finding of the indices is done correctly independently of the sorting of the input array elements.
//

int Functions::countArrayIndicesWithinBoundaries(RefArrayXd const array, double lowerBound, double upperBound)
{
    // At least one point is needed

    assert(array.size() >= 1);
    int binSize = 0;
    
    if (lowerBound < upperBound)
    {
        for (int i = 0; i < array.size(); ++i)
        {
            if ((array(i) >= lowerBound) && (array(i) <= upperBound))
                binSize++;
        }
    }

    return binSize;
}













// Functions::cubicSplineInterpolation()
// 
// PURPOSE:
//      This function computes a cubic spline interpolation of one-dimensional input data,
//      given a grid of new values we want the input data to be interpolated. The algorithm adopted is adapted
//      from Press W. et all, Numerical Recipes: The Art of Scientific Computing, 3rd Ed., 2007.
//
// INPUT:
//      observedAbscissa:                   an Eigen array containing the input data covariates to be used for the interpolation
//      observedOrdinate:                   an Eigen array containing the corresponding dependent variables of the input abscissa
//      interpolatedAbscissaUntruncated:    an Eigen array containing the values of the covariates for which we need to compute the
//                                          interpolated dependent variables. These values has not necessarily to be within the 
//                                          observed bounds
//
// OUTPUT:
//      An eigen array of doubles containing the newly computed ordinates (dependent variables) of the corresponding
//      input grid of covariates for which the interpolation was required.
//
// REMARKS:
//      It is required that all input abscissa arrays contain values sorted in ascending order.
//      No extrapolation can be done.
//

ArrayXd Functions::cubicSplineInterpolation(RefArrayXd const observedAbscissa, RefArrayXd const observedOrdinate, 
                                            RefArrayXd const interpolatedAbscissaUntruncated)
{
    // Number of data points
    
    int size = observedAbscissa.size();           
   

    // Number of interpolation grid points.
    
    int interpolatedSize = interpolatedAbscissaUntruncated.size();
    
   
    // Since the formula requires at least 2 data points, check if array size is not below 2,
    // if the interpolated grid has at least one point, and if input abscissa and ordinate 
    // have same number of elements.
    
    assert(size >= 2);
    assert(interpolatedSize >= 1);
    assert(observedOrdinate.size() == size);


    // Check if lower bound set by observed abscissa is respected by interpolated abscissa
   
    assert(observedAbscissa(0) <= interpolatedAbscissaUntruncated(0));
    

    // Compare upper bound of observed and interpolated abscissas and truncate the latter if it exceeds the observed upper limit

    double largestObservedAbscissa = observedAbscissa(size-1);
    double largestInterpolatedAbscissa = interpolatedAbscissaUntruncated(interpolatedSize-1);
    ArrayXd interpolatedAbscissa;

    if (largestObservedAbscissa < largestInterpolatedAbscissa)
    {
        // Since upper bound of observed abscissa is lower, and the routine is not doing any extrapolation, truncate the array of
        // interpolated abscissa at this upper bound.

        int extraSize = Functions::countArrayIndicesWithinBoundaries(interpolatedAbscissaUntruncated, largestObservedAbscissa, largestInterpolatedAbscissa);
        interpolatedSize = interpolatedSize - extraSize;
        interpolatedAbscissa = interpolatedAbscissaUntruncated.segment(0,interpolatedSize);
    }
    else
        interpolatedAbscissa = interpolatedAbscissaUntruncated;


    // Define some array differences

    ArrayXd differenceOrdinate = observedOrdinate.segment(1,size-1) - observedOrdinate.segment(0,size-1);
    ArrayXd differenceAbscissa = observedAbscissa.segment(1,size-1) - observedAbscissa.segment(0,size-1);
    ArrayXd differenceAbscissa2 = observedAbscissa.segment(2,size-2) - observedAbscissa.segment(0,size-2);


    // Lower bound condition for natural spline 

    vector<double> secondDerivatives(size);
    vector<double> partialSolution(size-1);
    secondDerivatives[0] = 0.0;
    partialSolution[0] = 0.0;
    
    
    // Do tridiagonal decomposition for computing second derivatives of observed ordinate
    
    ArrayXd sigma = differenceAbscissa.segment(0,size-2)/differenceAbscissa2.segment(0,size-2);
    double beta;

    
    // Forward computation of partial solutions in tridiagonal system

    for (int i = 1; i < size-1; ++i)
    {
        beta = sigma(i-1) * secondDerivatives[i-1] + 2.0;
        secondDerivatives[i] = (sigma(i-1) - 1.0)/beta;
        partialSolution[i] = differenceOrdinate(i)/differenceAbscissa(i) - differenceOrdinate(i-1)/differenceAbscissa(i-1);
        partialSolution[i] = (6.0*partialSolution[i]/differenceAbscissa2(i-1)-sigma(i-1)*partialSolution[i-1])/beta;
    }

    
    // Upper bound condition for natural spline
    
    secondDerivatives[size-1] = 0.0;


    // Backward substitution

    for (int k = (size-2); k >= 0; --k)
    {
        secondDerivatives[k] = secondDerivatives[k]*secondDerivatives[k+1]+partialSolution[k];
    }


    // Initialize arrays of differences in both ordinate and abscissa
    
    ArrayXd interpolatedOrdinate = ArrayXd::Zero(interpolatedSize);
    ArrayXd remainingInterpolatedAbscissa = interpolatedAbscissa;       // The remaining part of the array of interpolated abscissa
    int cumulatedBinSize = 0;                                           // The cumulated number of interpolated points from the beginning
    int i = 0;                                                          // Bin counter

    while ((i < size-1) && (cumulatedBinSize < interpolatedSize))
    {
        // Find which values of interpolatedAbscissa are containined within the selected bin of observedAbscissa.
        // Since elements in interpolatedAbscissa are monotonically increasing, we cut the input array each time
        // we identify the elements of the current bin. This allows to speed up the process.

        double lowerAbscissa = observedAbscissa(i);
        double upperAbscissa = observedAbscissa(i+1);


        // Find total number of interpolated points falling in the current bin

        int binSize = Functions::countArrayIndicesWithinBoundaries(remainingInterpolatedAbscissa, lowerAbscissa, upperAbscissa);
      

        // Do interpolation only if interpolated points are found within the bin

        if (binSize > 0)
        {
            double lowerOrdinate = observedOrdinate(i);
            double upperOrdinate = observedOrdinate(i+1);
            double denominator = differenceAbscissa(i);
            double upperSecondDerivative = secondDerivatives[i+1];
            double lowerSecondDerivative = secondDerivatives[i];
            ArrayXd interpolatedAbscissaInCurrentBin = remainingInterpolatedAbscissa.segment(0, binSize);
            ArrayXd interpolatedOrdinateInCurrentBin = ArrayXd::Zero(binSize);


            // Compute coefficients for cubic spline interpolation function

            ArrayXd a = (upperAbscissa - interpolatedAbscissaInCurrentBin) / denominator;
            ArrayXd b = 1.0 - a;
            ArrayXd c = (1.0/6.0) * (a.cube() - a)*denominator*denominator;
            ArrayXd d = (1.0/6.0) * (b.cube() - b)*denominator*denominator;
            interpolatedOrdinateInCurrentBin = a*lowerOrdinate + b*upperOrdinate + c*lowerSecondDerivative + d*upperSecondDerivative;

                
            // Merge bin ordinate into total array of ordinate
        
            interpolatedOrdinate.segment(cumulatedBinSize, binSize) = interpolatedOrdinateInCurrentBin;


            // Reduce size of array remainingInterpolatedAbscissa by binSize elements and initialize the array
            // with remaining part of interpolatedAbscissa
        
            int currentRemainingSize = interpolatedSize - cumulatedBinSize;
            remainingInterpolatedAbscissa.resize(currentRemainingSize - binSize);
            cumulatedBinSize += binSize;
            remainingInterpolatedAbscissa = interpolatedAbscissa.segment(cumulatedBinSize, interpolatedSize-cumulatedBinSize);
        }  


        // Move to next bin

        ++i;
    }

    return interpolatedOrdinate;
}
