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
     
        virtual double logDensity(RefArrayXd x, const bool includeConstantTerm=false);
        virtual void draw(RefArrayXXd sample);
        virtual void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood);
        virtual bool isUniformPrior();


    private:
        
        vector<normal_distribution<>> normalDistributionVector;
        normal_distribution<> normal;
        ArrayXd mean;
        ArrayXd standardDeviation;


}; // END class NormalPrior

#endif
