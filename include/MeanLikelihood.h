// Derived class for mean likelihood computations
// Created by Enrico Corsaro @ IvS - 19 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "MeanLikelihood.h"
// Implementations contained in "MeanLikelihood.cpp"


#ifndef MEANLIKELIHOOD_H
#define MEANLIKELIHOOD_H

#include <cmath>
#include <iostream>
#include "Likelihood.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class MeanLikelihood : public Likelihood
{

    public:

        MeanLikelihood(const RefArrayXd covariates, const RefArrayXd observations, const RefArrayXd uncertainties, Model &model);
        ~MeanLikelihood();
        ArrayXd getNormalizedUncertainties();

        virtual double logValue(RefArrayXd modelParameters);

    private:
        
        ArrayXd mormalizedUncertainties;
}; // END class MeanLikelihood

#endif
