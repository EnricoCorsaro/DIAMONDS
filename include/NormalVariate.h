
#ifndef NORMALVARIATE_H
#define NORMALVARIATE_H

#include <vector>
#include <random>
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