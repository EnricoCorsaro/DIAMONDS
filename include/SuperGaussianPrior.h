// Derived class for Super-Gaussian prior computations
// Created by Enrico Corsaro @ IvS - 6 August 2013
// e-mail: emncorsaro@gmail.com
// Header file "SuperGaussianPrior.h"
// Implementations contained in "SuperGaussianPrior.cpp"


#ifndef SUPERGAUSSIANPRIOR_H
#define SUPERGAUSSIANPRIOR_H

#include <iostream>
#include "Prior.h"


using namespace std;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class SuperGaussianPrior : public Prior
{
    public:

        SuperGaussianPrior(const RefArrayXd center, const RefArrayXd sigma, const RefArrayXd widthOfPlateau);
        ~SuperGaussianPrior();

        ArrayXd getCenter();
        ArrayXd getSigma();
        ArrayXd getWidthOfPlateau();
     
        virtual double logDensity(RefArrayXd const x, const bool includeConstantTerm = false);
        virtual bool drawnPointIsAccepted(RefArrayXd const drawnPoint);
        virtual void draw(RefArrayXXd drawnSample);
        virtual void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood);
        virtual void writeHyperParametersToFile(string fullPath);


    private:
        
        vector<normal_distribution<>> normalDistributionVector;
        normal_distribution<> normal;
        uniform_real_distribution<> uniform;
        ArrayXd center;
        ArrayXd sigma;
        ArrayXd widthOfPlateau;
        ArrayXd halfWidthOfPlateau;
        ArrayXd plateauArea;
        ArrayXd tailsArea;
        ArrayXd totalArea;
        ArrayXXd normalizedAreas;

}; // END class SuperGaussianPrior

#endif
