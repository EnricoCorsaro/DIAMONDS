// Derived class for exponential likelihood computations
// Created by Enrico Corsaro @ IvS - 5 June 2013
// e-mail: emncorsro@gmail.com
// Header file "ExponentialLikelihood.h"
// Implementations contained in "ExponentialLikelihood.cpp"


#ifndef EXPONENTIALLIKELIHOOD_H
#define EXPONENTIALLIKELIHOOD_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "Likelihood.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class ExponentialLikelihood : public Likelihood
{

    public:

        ExponentialLikelihood(const RefArrayXd observations, Model &model);
        ~ExponentialLikelihood();

        virtual double logValue(RefArrayXd const modelParameters);


    private:

}; // END class ExponentialLikelihood

#endif
