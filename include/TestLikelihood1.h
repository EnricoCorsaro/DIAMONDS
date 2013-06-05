// Class for Toy Likelihood example #1 by Feroz et al. 2008, MNRAS, 384, 449
// Created by Enrico Corsaro @ IvS - 28 May 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "TestLikelihood1.h"
// Implementations contained in "TestLikelihood1.cpp"


#ifndef TESTLIKELIHOOD1_H
#define TESTLIKELIHOOD1_H

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


class TestLikelihood1 : public Likelihood
{

    public:

        TestLikelihood1(const RefArrayXd observations, Model &model);
        ~TestLikelihood1();

        virtual double logValue(RefArrayXd nestedSampleOfParameters);


    private:

}; // END class TestLikelihood1

#endif
