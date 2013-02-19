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
#include <sstream>
#include <Eigen/Core>


using namespace std;

typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
using Eigen::ArrayXd;


namespace MathExtra
{

    const double PI = 4.0 * atan(1.0);

    // Functional Profiles
    void lorentzProfile(RefArrayXd y, const RefArrayXd x, double x0, double gamma, double amp = 1);
    double logGaussProfile(double x, double mu, double sigma, double amp = 1);
    void logGaussProfile(RefArrayXd y, const RefArrayXd x, const double mu, const double sigma, const double amp = 1);
    
    // Likelihood functions
    double logGaussLikelihood(const RefArrayXd observations, const RefArrayXd predictions, const RefArrayXd uncertainties);
    
    // Array functions
    inline double product(const vector<double> &vec);
    inline double sum(const vector<double> &vec);
    double logExpSum(double x, double y);

} // END namespace MathExtra
#endif
