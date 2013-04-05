// Derived class for normal prior computations
// Created by Enrico Corsaro @ IvS - 4 April 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "NormalPrior.h"
// Implementations contained in "NormalPrior.cpp"


#ifndef NORMALPRIOR_H
#define NORMALPRIOR_H

#include <iostream>
#include "Prior.h"


using namespace std;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class NormalPrior : public Prior
{
    public:

        NormalPrior(const RefArrayXd mean, const RefArrayXd standardDeviation);
        ~NormalPrior();

        ArrayXd getMean();
        ArrayXd getStandardDeviation();
        
        virtual double getNormalizingFactor();
        virtual void draw(RefArrayXXd nestedSampleOfParameters, const int Nobjects);
        virtual void drawWithConstraint(RefArrayXd nestedSampleOfParameters, Likelihood &likelihood);

    private:
        
        vector<normal_distribution<>> normalDistributionVector;
        normal_distribution<> normal;
        mt19937 engine;
        ArrayXd mean;
        ArrayXd standardDeviation;

}; // END class NormalPrior

#endif
