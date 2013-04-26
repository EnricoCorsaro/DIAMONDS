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

        virtual double getNormalizingFactor();
        virtual void draw(RefArrayXXd nestedSampleOfParameters, const int Nobjects);
        virtual void drawWithConstraint(RefArrayXd nestedSampleOfParameters, Likelihood &likelihood);
        virtual bool pointIsRejected(RefArrayXXd drawnSampleOfParameters);


    private:

        uniform_real_distribution<> uniform;
        ArrayXd minima;
        ArrayXd maxima;


}; // END class UniformPrior

#endif
