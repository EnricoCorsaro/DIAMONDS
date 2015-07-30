// Derived class for uniform prior computations
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 15 February 2013
// e-mail: emncorsaro@gmail.com
// Header file "UniformPrior.h"
// Implementations contained in "UniformPrior.cpp"


#ifndef UNIFORMPRIOR_H
#define UNIFORMPRIOR_H

#include <iostream>
#include <limits>
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

        virtual double logDensity(RefArrayXd const x, const bool includeConstantTerm = false);
        virtual bool drawnPointIsAccepted(RefArrayXd const drawnPoint);
        virtual void draw(RefArrayXXd drawnSample);
        virtual void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood);
        virtual void writeHyperParametersToFile(string fullPath);


    private:

        uniform_real_distribution<> uniform;
        ArrayXd minima;
        ArrayXd maxima;


}; // END class UniformPrior

#endif
