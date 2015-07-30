// Derived class for zero prior objects (used only in the merger)
// This class returns zero (or false) values for any input parameters
// Created by Enrico Corsaro @ IvS - 6 June 2014
// e-mail: emncorsaro@gmail.com
// Header file "ZeroPrior.h"
// Implementations contained in "ZeroPrior.cpp"


#ifndef ZEROPRIOR_H
#define ZEROPRIOR_H

#include <iostream>
#include <limits>
#include "Prior.h"


using namespace std;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class ZeroPrior : public Prior
{

    public:

        ZeroPrior(const int Ndimensions);
        ~ZeroPrior();

        virtual double logDensity(RefArrayXd const x, const bool includeConstantTerm);
        virtual bool drawnPointIsAccepted(RefArrayXd const drawnPoint);
        virtual void draw(RefArrayXXd drawnSample){};
        virtual void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood){};
        virtual void writeHyperParametersToFile(string fullPath){};


    private:


}; // END class ZeroPrior

#endif
