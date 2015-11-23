// Derived class for uniform prior computations
// Created by Enrico Corsaro - November 2015
// e-mail: emncorsaro@gmail.com
// Header file "GridUniformPrior.h"
// Implementations contained in "GridUniformPrior.cpp"


#ifndef GRIDUNIFORMPRIOR_H
#define GRIDUNIFORMPRIOR_H

#include <iostream>
#include <limits>
#include "Prior.h"


using namespace std;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class GridUniformPrior : public Prior
{

    public:

        GridUniformPrior(const RefArrayXd startingCoordinate, const RefArrayXd NgridPoints, const RefArrayXd separation, 
                         const RefArrayXd tolerance);
        ~GridUniformPrior();

        ArrayXd getStartingCoordinate();
        ArrayXd getNgridPoints();
        ArrayXd getSeparation();
        ArrayXd getTolerance();

        virtual double logDensity(RefArrayXd const x, const bool includeConstantTerm = false);
        virtual bool drawnPointIsAccepted(RefArrayXd const drawnPoint);
        virtual void draw(RefArrayXXd drawnSample);
        virtual void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood);
        virtual void writeHyperParametersToFile(string fullPath);


    private:

        uniform_real_distribution<> uniform;
        vector<uniform_int_distribution<>> uniformIntegerVector;
        ArrayXd startingCoordinate;
        ArrayXd NgridPoints;
        ArrayXd separation;
        ArrayXd tolerance;

};

#endif
