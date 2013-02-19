// Derived class for normal likelihood computations
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 15 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "NormalLikelihood.h"
// Implementations contained in "NormalLikelihood.cpp"


#ifndef NORMALLIKELIHOOD_H
#define NORMALLIKELIHOOD_H

#include <cmath>
#include <iostream>
#include "Likelihood.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class NormalLikelihood : public Likelihood
{

    public:

        NormalLikelihood(const RefArrayXd covariates, const RefArrayXd observations, const RefArrayXd uncertainties, Model &model);
        ~NormalLikelihood();
        ArrayXd getCovariates();
        ArrayXd getObservations();
        ArrayXd getUncertainties();

        virtual double logValue(RefArrayXd modelParameters);

    private:

}; // END class NormalLikelihood

#endif
