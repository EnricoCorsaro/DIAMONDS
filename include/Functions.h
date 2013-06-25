// Class for implementing some extra functions useful 
// computations of models and likelihoods.
// Created by Enrico Corsaro @ IvS - 22 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
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
    double logGaussProfile(const double covariate, const double mu = 0.0, 
                           const double sigma = 1.0, const double amplitude = 1.0);
    void logGaussProfile(RefArrayXd predictions, const RefArrayXd covariates, const double mu = 0.0, 
                         const double sigma = 1.0, const double amplitude = 1.0);
    
    
    // Likelihood functions
    
    double logGaussLikelihood(const RefArrayXd observations, const RefArrayXd predictions, const RefArrayXd uncertainties);
    
    
    // Matrix algebra functions

    void clusterCovariance(const RefArrayXXd clusterSample, RefArrayXXd covarianceMatrix, 
                           RefArrayXd centerCoordinates);
    void selfAdjointMatrixDecomposition(const RefArrayXXd covarianceMatrix, RefArrayXd eigenvalues, 
                                        RefArrayXXd eigenvectorsMatrix);
    

    // Array algebra functions
    
    inline double product(const vector<double> &vec);
    inline double sum(const vector<double> &vec);
    double logExpSum(const double x, const double y);
    double logExpDifference(const double x, const double y);
    void sortElementsDouble(RefArrayXd array1, RefArrayXd array2);   // TODO: replace by much more efficient algorithms 
    void sortElementsInt(vector<int> &array1, RefArrayXd array2);    // Idem


    // Utility functions

    template <typename Type>
    vector<int> argsort(const vector<Type> &myVector);


} // END namespace Functions





// Functions::argsort()
//
// PURPOSE: 
//      Mimicks Numpy's argsort() function. Given a vector of any type, it returns a
//      vector of integer indices such that
//        myVector[indices[0]] <= myVector[indices[1]] <= ... <= myVector[indices[N-1]]
//
// INPUT:
//      myVector[0..N-1]: vector which you want to argsort
//
// OUTPUT:
//      indices[0-N-1]: vector with integer indices
//
//
// REMARKS:
//      Because of the template, this function needs to be in the header file
// 

template <typename Type>
vector<int> Functions::argsort(const vector<Type> &myVector)
{
    vector<int> indices(myVector.size());

    // Generate the indices from 0 to N-1

    iota(begin(indices), end(indices), 0);

    // Use the sort() function of the STL, but with a lambda function
    // that compares the elements of myVector

    sort(begin(indices), end(indices), [&myVector] (int i, int j) {return myVector[i] < myVector[j];} );

    return indices;
}




#endif
