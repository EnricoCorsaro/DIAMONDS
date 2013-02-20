// Derived class for uniform prior computations
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 15 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "UniformPrior.h"
// Implementations contained in "UniformPrior.cpp"


#ifndef UNIFORMPRIOR_H
#define UNIFORMPRIOR_H

#include <iostream>
#include "Prior.h"


using namespace std;
using Eigen:ArrayXd;
using Eigen:ArrayXXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class UniformPrior : public Prior
{

    public:

        UniformPrior(const RefArrayXXd boundaries, const int Nobjects);
        ~UniformPrior();

        ArrayXXd getBoundaries();
        double getUniformFactor();

        virtual void draw(RefArrayXXd nestedParameters);
        virtual void drawWithConstraint(RefArrayXd nestedParameters, Likelihood &likelihood);

    private:

        uniform_real_distribution<> uniform;
        mt19937 engine;
        ArrayXXd boundaries;
        double uniformFactor;

}; // END class UniformPrior

#endif
