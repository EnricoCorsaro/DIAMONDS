// Class for implementing object construction from prior and likelihood distributions.
// This class is to be used within the NestedSampler class.
// Created by Joris De Ridder @ IvS - 01 February 2013
// e-mail: joris.deridder@ster.kuleuven.be
// Header file "NormalVariate.h"


#ifndef NORMALVARIATE_H
#define NORMALVARIATE_H

#include <vector>
#include <random>
#include <ctime>
#include "RandomVariate.h"
#include "MathExtra.h"


using namespace std;

class NormalVariate : public RandomVariate
{
    public:
    
        NormalVariate(double mu, double sigma);
        ~NormalVariate();
        virtual void drawNestedValues(vector<double> &values, vector<double> &logDensities, int Nvalues);
        virtual void drawNestedValueWithConstraint(double &value, double &logDensity, double logDensityConstraint);
        
    protected:
    
    private:
    
        uniform_real_distribution<> uniform;
        mt19937 engine;
        double mu;
        double sigma;
    
};


#endif
