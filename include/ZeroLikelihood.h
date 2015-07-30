// Derived class for zero likelihood objects
// This likelihood returns zero for any input set of observations and model parameters.
// Created by Enrico Corsaro @ IvS - 9 May 2014
// e-mail: emncorsaro@gmail.com
// Header file "ZeroLikelihood.h"
// Implementations contained in "ZeroLikelihood.cpp"


#ifndef ZEROLIKELIHOOD_H
#define ZEROLIKELIHOOD_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "Likelihood.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class ZeroLikelihood : public Likelihood
{

    public:

        ZeroLikelihood(const RefArrayXd observations, Model &model);
        ~ZeroLikelihood();

        virtual double logValue(RefArrayXd const modelParameters) override;


    private:

};

#endif
