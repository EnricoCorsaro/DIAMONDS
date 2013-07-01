// Abstract base class for prior computations.
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 15 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "Prior.h"
// Implementations contained in "Prior.cpp"


#ifndef PRIOR_H
#define PRIOR_H

#include <random>
#include <ctime>
#include <vector>
#include <cstdlib>
#include <Eigen/Core>
#include "Likelihood.h"
#include "Functions.h"

using namespace std;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class Prior
{

    public:

        Prior(const int Ndimensions);
        ~Prior();
        int getNdimensions();
        
        virtual double logDensity(RefArrayXd x, const bool includeConstantTerm=false) = 0;
        virtual void draw(RefArrayXXd sample) = 0;
        virtual void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood) = 0;

        const double minusInfinity;

    protected:
        
        int Ndimensions;
        mt19937 engine;

    
    private:
    

}; // END class Prior

#endif
