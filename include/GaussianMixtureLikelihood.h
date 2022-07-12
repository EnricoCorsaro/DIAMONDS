// Derived class for a Gaussian Mixture Likelihood
// Created by Enrico Corsaro @ OACT - September 2020
// e-mail: emncorsaro@gmail.com
// Header file "GaussianMixtureLikelihood.h"
// Implementations contained in "GaussianMixtureLikelihood.cpp"


#ifndef GAUSSIANMIXTURELIKELIHOOD_H
#define GAUSSIANMIXTURELIKELIHOOD_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "Likelihood.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class GaussianMixtureLikelihood : public Likelihood
{

    public:

        GaussianMixtureLikelihood(const RefArrayXd observations, const RefArrayXd uncertainties, Model &model);
        ~GaussianMixtureLikelihood();
        ArrayXd getUncertainties();

        virtual double logValue(RefArrayXd const modelParameters);


    private:

        ArrayXd uncertainties;
        int Npoints;

}; 

#endif
