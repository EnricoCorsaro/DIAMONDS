// Class for Toy Likelihood example #2 by Feroz et al. 2009, MNRAS, 398, 1601
// Created by Enrico Corsaro @ IvS - 29 May 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "TestLikelihood2.h"
// Implementations contained in "TestLikelihood2.cpp"


#ifndef TESTLIKELIHOOD2_H
#define TESTLIKELIHOOD2_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "Functions.h"
#include <Eigen/Core>
#include "Likelihood.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class TestLikelihood2 : public Likelihood
{

    public:

        TestLikelihood2(const RefArrayXd observations, Model &model);
        ~TestLikelihood2();

        virtual double logValue(RefArrayXd modelParameters);


    private:

};

#endif
