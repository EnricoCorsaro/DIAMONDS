// Derived class for uniform prior computations
// Created by Enrico Corsaro @ IvS - 28 October 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "DeltaPrior.h"
// Implementations contained in "DeltaPrior.cpp"


#ifndef DELTAPRIOR_H
#define DELTAPRIOR_H

#include <iostream>
#include <limits>
#include "Prior.h"


using namespace std;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class DeltaPrior : public Prior
{

    public:

        DeltaPrior(const RefArrayXd constantParameters);
        ~DeltaPrior();

        ArrayXd getConstants();

        virtual double logDensity(RefArrayXd x, const bool includeConstantTerm=false);
        virtual void draw(RefArrayXXd drawnSample);


    private:

        ArrayXd constantParameters;

};

#endif
