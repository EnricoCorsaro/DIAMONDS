// Derived class for normal likelihood computations
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 15 February 2013
// e-mail: emncorsaro@gmail.com
// Header file "NormalLikelihood.h"
// Implementations contained in "NormalLikelihood.cpp"


#ifndef NORMALLIKELIHOOD_H
#define NORMALLIKELIHOOD_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "Likelihood.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class NormalLikelihood : public Likelihood
{

    public:

        NormalLikelihood(const RefArrayXd observations, const RefArrayXd uncertainties, Model &model);
        ~NormalLikelihood();
        ArrayXd getUncertainties();

        virtual double logValue(RefArrayXd const modelParameters);


    private:

        ArrayXd uncertainties;

}; 

#endif
