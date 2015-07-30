// Derived class for uniform prior computations
// Created by Enrico Corsaro @ IvS - 6 March 2013
// e-mail: emncorsaro@gmail.com
// Header file "GridPrior.h"
// Implementations contained in "GridPrior.cpp"


#ifndef GRIDPRIOR_H
#define GRIDPRIOR_H

#include <iostream>
#include <limits>
#include "Prior.h"


using namespace std;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class GridPrior : public Prior
{

    public:

        GridPrior(const RefArrayXd width, const RefArrayXd separation, const RefArrayXd startingCoordinate, const RefArrayXd Nsteps);
        ~GridPrior();

        ArrayXd getWidth();
        ArrayXd getSeparation();
        ArrayXd getStartingCoordinate();
        ArrayXd getNsteps();

        virtual double logDensity(RefArrayXd const x, const bool includeConstantTerm = false);
        virtual bool drawnPointIsAccepted(RefArrayXd const drawnPoint);
        virtual void draw(RefArrayXXd drawnSample);
        virtual void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood);
        virtual void writeHyperParametersToFile(string fullPath);


    private:

        uniform_real_distribution<> uniform;
        vector<uniform_int_distribution<>> uniformIntegerVector;
        ArrayXd width;
        ArrayXd separation;
        ArrayXd startingCoordinate;
        ArrayXd Nsteps;

};

#endif
