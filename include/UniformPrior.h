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
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class UniformPrior : public Prior
{

    public:

        UniformPrior(const RefArrayXd minima, const RefArrayXd maxima);
        ~UniformPrior();

        ArrayXd getMinima();
        ArrayXd getMaxima();
        double getUniformFactor();

        virtual void draw(RefArrayXXd nestedParameters, const int Nobjects);
        virtual void drawWithConstraint(RefArrayXd nestedParameters, Likelihood &likelihood);

    private:

        uniform_real_distribution<> uniform;
        mt19937 engine;
        ArrayXd minima;
        ArrayXd maxima;
        double uniformFactor;

}; // END class UniformPrior

#endif
