// Abstract base class for prior computations.
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 15 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "Prior.h"
// Implementations contained in "Prior.cpp"


#ifndef PRIOR_H
#define PRIOR_H

#include <random>
#include <cstdlib>
#include <Eigen/Core>
#include "Likelihood.h"

using namespace std;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class Prior
{

    public:

        Prior(const int Ndimensions);
        ~Prior();
        int getNdimensions();


        // Pure virtual functions implemented in derived classes
        
        virtual void draw(RefArrayXXd nestedParameters, const int Nobjects) = 0;
        virtual void drawWithConstraint(RefArrayXd nestedParameters, Likelihood &likelihood) = 0;

    protected:
        
        int Ndimensions;

    private:
    
}; // END class Prior

#endif
