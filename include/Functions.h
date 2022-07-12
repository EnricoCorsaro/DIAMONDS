// Class for implementing some extra functions useful 
// computations of models and likelihoods.
// Created by Enrico Corsaro @ IvS - 22 January 2013
// e-mail: emncorsaro@gmail.com
// Header file "Functions.h"
// Implementation contained in "Functions.cpp"


#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <cassert>
#include <numeric>
#include <functional>
#include <random>
#include <ctime>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <Eigen/Dense>
#include "Metric.h"

#define SWAPDOUBLE(a,b) {double copy; copy = a, a = b, b = copy;}
#define SWAPINT(a,b) {int copy; copy = a, a = b, b = copy;}


using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXi> RefArrayXi;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;

namespace Functions
{

    const double PI = 4.0 * atan(1.0);


    // Profile functions
    
    void lorentzProfile(RefArrayXd predictions, const RefArrayXd covariates, 
                        const double centroid = 0.0, const double amplitude = 1.0, const double gamma = 1.0);
    void modeProfile(RefArrayXd predictions, const RefArrayXd covariates, 
                     const double centroid = 0.0, const double height = 1.0, const double linewidth = 1.0);
    void modeProfileWithAmplitude(RefArrayXd predictions, const RefArrayXd covariates, 
                                  const double centroid = 0.0, const double amplitude = 1.0, const double linewidth = 1.0);
    void modeProfileAsymmetricWithAmplitude(RefArrayXd predictions, const RefArrayXd covariates, const double centroid = 0.0, 
                                            const double amplitude = 1.0, const double linewidth = 1.0, const double asymmetry = 0.0);
    void modeProfileSinc(RefArrayXd predictions, const RefArrayXd covariates, 
                     const double centroid = 0.0, const double height = 1.0, const double resolution = 1.0);
    double logGaussProfile(const double covariate, const double mu = 0.0, 
                           const double sigma = 1.0, const double amplitude = 1.0);
    void logGaussProfile(RefArrayXd predictions, const RefArrayXd covariates, const double mu = 0.0, 
                         const double sigma = 1.0, const double amplitude = 1.0);
    
    
    // Likelihood functions
    
    double logGaussLikelihood(const RefArrayXd observations, const RefArrayXd predictions, const RefArrayXd uncertainties);
    
    
    // Matrix algebra functions

    void clusterCovariance(RefArrayXXd const clusterSample, RefArrayXXd covarianceMatrix, 
                           RefArrayXd centerCoordinates);
    bool selfAdjointMatrixDecomposition(RefArrayXXd const covarianceMatrix, RefArrayXd eigenvalues, 
                                        RefArrayXXd eigenvectorsMatrix);
    

    // Array manipulation functions
    
    inline double product(const vector<double> &vec);
    inline double sum(const vector<double> &vec);
    double logExpSum(const double x, const double y);
    double logExpDifference(const double x, const double y);
    void topDownMerge(RefArrayXd array1, RefArrayXd arrayCopy1,         // Only used within topDownMergeSort
                      RefArrayXd array2, RefArrayXd arrayCopy2, 
                      int beginIndex, int middleIndex, int endIndex);
    void topDownSplitMerge(RefArrayXd array1, RefArrayXd arrayCopy1,    // Only used within topDownMergeSort
                           RefArrayXd array2, RefArrayXd arrayCopy2, 
                           int beginIndex, int endIndex);
    void topDownMergeSort(RefArrayXd array1, RefArrayXd array2);        // Mergesort algorithm to sort array1 in ascending order and array2 according to array1
    void sortElementsDouble(RefArrayXd array1, RefArrayXd array2);   
    void sortElementsInt(vector<int> &array1, RefArrayXd array2);


    vector<int> findArrayIndicesWithinBoundaries(RefArrayXd const array, double lowerBound, double upperBound);
    int countArrayIndicesWithinBoundaries(RefArrayXd const array, double lowerBound, double upperBound);
    ArrayXd cubicSplineInterpolation(RefArrayXd const observedAbscissa, RefArrayXd const observedOrdinate, 
                                     RefArrayXd const interpolatedAbscissaUntruncated);


    // Utility functions

    long double factorial(int number);
    template <typename Type>
    vector<int> argsort(const vector<Type> &myVector);

}









// Functions::argsort()
//
// PURPOSE: 
//      Mimicks Numpy's argsort() function. Given a vector of any type, it returns a
//      vector of integer indices such that:
//      myVector[indices[0]] <= myVector[indices[1]] <= ... <= myVector[indices[N-1]]
//
// INPUT:
//      myVector[0..N-1]: vector which you want to sort in ascending order
//
// OUTPUT:
//      indices[0-N-1]: vector with integer indices corresponding to sorted elements
//                      of the input vector
//
// REMARKS:
//      Because of the template, this function needs to be in the header file.
// 

template <typename Type>
vector<int> Functions::argsort(const vector<Type> &myVector)
{
    // Create a vector of integers to contain the indices of the sorted elements
    
    vector<int> indices(myVector.size());


    // Generate the indices from 0 to myVector.size()-1
    // std::iota is implemented in C++11 standard. 
    // It requires the starting and ending indices of the array to fill with increasing integers
    // and the starting value, which is incremented at each subsequent element (++value).

    iota(begin(indices), end(indices), 0);


    // Use the sort() function implemented in algorithm library, but with a lambda function
    // that compares the elements of myVector.

    sort(begin(indices), end(indices), [&myVector] (int i, int j) {return myVector[i] < myVector[j];} );

    return indices;
}


#endif
