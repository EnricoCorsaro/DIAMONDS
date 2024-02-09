// Derived class for normal likelihood computations of a generalized model
// Created by Enrico Corsaro @ OACT - December 2022
// e-mail: enrico.corsaro@inaf.it
// Header file "GeneralizedNormalLikelihood.h"
// Implementations contained in "GeneralizedNormalLikelihood.cpp"


#ifndef GENERALIZEDNORMALLIKELIHOOD_H
#define GENERALIZEDNORMALLIKELIHOOD_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "Likelihood.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class GeneralizedNormalLikelihood : public Likelihood
{

    public:

        GeneralizedNormalLikelihood(const RefArrayXd observations, const RefArrayXd observationsUncertainty, 
                                    Model &model);
        ~GeneralizedNormalLikelihood();
        ArrayXd getObservationsUncertainty();

        virtual double logValue(RefArrayXd const modelParameters);


    private:

        ArrayXd observationsUncertainty;

}; 

#endif
