// Derived class for normal likelihood computations of a multilinear model
// Created by Enrico Corsaro @ OACT - December 2017
// e-mail: emncorsaro@gmail.com
// Header file "MultiLinearNormalLikelihood.h"
// Implementations contained in "MultiLinearNormalLikelihood.cpp"


#ifndef MULTILINEARNORMALLIKELIHOOD_H
#define MULTILINEARNORMALLIKELIHOOD_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "Likelihood.h"
#include "MultiLinearModel.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class MultiLinearNormalLikelihood : public Likelihood
{

    public:

        MultiLinearNormalLikelihood(const RefArrayXd observations, const RefArrayXd covariatesUncertainties, 
                                    const RefArrayXd observationsUncertainty, MultiLinearModel &model);
        ~MultiLinearNormalLikelihood();
        ArrayXd getCovariatesUncertainties();
        ArrayXd getObservationsUncertainty();

        virtual double logValue(RefArrayXd const modelParameters);


    private:

        ArrayXd covariatesUncertainties;
        ArrayXd observationsUncertainty;
        int Nobservables;
        int Npoints;

}; 

#endif
