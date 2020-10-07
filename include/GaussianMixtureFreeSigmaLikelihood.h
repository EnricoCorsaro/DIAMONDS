// Derived class for normal likelihood computations
// Created by Enrico Corsaro @ OACT - September 2020
// e-mail: emncorsaro@gmail.com
// Header file "GaussianMixtureFreeSigmaLikelihood.h"
// Implementations contained in "GaussianMixtureFreeSigmaLikelihood.cpp"


#ifndef GAUSSIANMIXTUREFREESIGMALIKELIHOOD_H
#define GAUSSIANMIXTUREFREESIGMALIKELIHOOD_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "Likelihood.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class GaussianMixtureFreeSigmaLikelihood : public Likelihood
{

    public:

        GaussianMixtureFreeSigmaLikelihood(const RefArrayXd observations, Model &model);
        ~GaussianMixtureFreeSigmaLikelihood();

        virtual double logValue(RefArrayXd const modelParameters);


    private:

}; 

#endif
