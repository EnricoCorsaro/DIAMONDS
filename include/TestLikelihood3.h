// Class for Toy Likelihood example #2 by Feroz et al. 2008, MNRAS, 384, 449
// Created by Enrico Corsaro @ IvS - 30 May 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "TestLikelihood3.h"
// Implementations contained in "TestLikelihood3.cpp"


#ifndef TESTLIKELIHOOD3_H
#define TESTLIKELIHOOD3_H

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


class TestLikelihood3 : public Likelihood
{

    public:

        TestLikelihood3(const RefArrayXd observations, Model &model);
        ~TestLikelihood3();

        virtual double logValue(RefArrayXd nestedSampleOfParameters);


    private:

}; // END class TestLikelihood3

#endif
