// Class for implementing some extra functions useful for Peak Bagging fitting process
// Created by Enrico Corsaro @ IvS - 22 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "MathExtra.h"


#ifndef MATHEXTRA_H
#define MATHEXTRA_H


#include <cmath>
#include <numeric>
#include <functional>
#include <vector>
#include <iostream>
#include <sstream>


const double PI = 4.0 * atan(1.0);

using namespace std;

class MathExtra
{
    public:
    
        vector<double> y;
        void lorentzProfile(vector<double> &x, double x0, double gamma, double amp = 1);
        void gaussProfile(vector<double> &x, double x0, double sigma, double amp = 1);
        double gaussLikelihood(vector<double> &x_obs, vector<double> &x_theor, vector<double> &sigma);
        double logGaussLikelihood(vector<double> &x_obs, vector<double> &x_theor, vector<double> &sigma);

        inline double product(const vector<double> &vec);
        inline double sum(const vector<double> &vec);
        inline double logExpSum(const double x, const double y);
    
    protected:
    
    private:
    
};


#endif
