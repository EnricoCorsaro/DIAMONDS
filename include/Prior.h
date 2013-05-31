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

        Prior(const int Ndimensions, const bool uniformFlag);
        ~Prior();
        int getNdimensions();
        bool priorIsUniform();     // true if uniform, false otherwise
        
        virtual double getNormalizingFactor() = 0;
        virtual void draw(RefArrayXXd nestedSampleOfParameters) = 0;
        virtual void drawWithConstraint(RefArrayXd nestedSampleOfParameters, Likelihood &likelihood) = 0;
        virtual bool pointIsRejected(RefArrayXXd drawnSampleOfParameters) = 0;


    protected:
        
        bool uniformFlag;
        int Ndimensions;
        double normalizingFactor;
        mt19937 engine;

    
    private:
    

}; // END class Prior

#endif
