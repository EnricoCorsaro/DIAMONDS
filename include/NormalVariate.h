// Class for implementing object construction from prior and likelihood distributions.
// This class is to be used within the NestedSampler class.
// Created by Joris De Ridder @ IvS - 01 February 2013
// e-mail: joris.deridder@ster.kuleuven.be
// Header file "NormalVariate.h"


#ifndef NORMALVARIATE_H
#define NORMALVARIATE_H

#include <random>
#include <ctime>
#include <Eigen/Core>
#include "RandomVariate.h"
#include "MathExtra.h"


using namespace std;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class NormalVariate : public RandomVariate
{
    public:
    
        NormalVariate(double mu, double sigma);
        ~NormalVariate();
        virtual void drawNestedValues(RefArrayXXd values, RefArrayXd logDensities, int Nvalues);
        virtual void drawNestedValueWithConstraint(RefArrayXd value, double &logDensity, double logDensityConstraint);
                
    protected:
    
    private:
    
        uniform_real_distribution<> uniform;
        mt19937 engine;
        double mu;
        double sigma;
    
};


#endif
