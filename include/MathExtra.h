// Class for implementing some extra functions useful 
// computations of models and likelihoods.
// Created by Enrico Corsaro @ IvS - 22 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "MathExtra.h"
// Implementation contained in "MathExtra.cpp"


#ifndef MATHEXTRA_H
#define MATHEXTRA_H

#include <cmath>
#include <numeric>
#include <functional>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <Eigen/Core>

#define SWAP(a,b) {double copy; copy = a, a = b, b = copy;}


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


namespace MathExtra
{

    const double PI = 4.0 * atan(1.0);

    // Functional Profiles
    void lorentzProfile(RefArrayXd predictions, const RefArrayXd covariates, const double centroid = 0.0, const double amplitude = 1.0, const double gamma = 1.0);
    double logGaussProfile(const double covariate, const double mu = 0.0, const double sigma = 1.0, const double amplitude = 1.0);
    void logGaussProfile(RefArrayXd y, const RefArrayXd x, const double mu = 0.0, const double sigma = 1.0, const double amplitude = 1.0);
    
    // Likelihood functions
    double logGaussLikelihood(const RefArrayXd observations, const RefArrayXd predictions, const RefArrayXd uncertainties);
    
    // Array functions
    inline double product(const vector<double> &vec);
    inline double sum(const vector<double> &vec);
    double logExpSum(double x, double y);
    void sortElements(RefArrayXd array1, RefArrayXd array2);

} // END namespace MathExtra
#endif
