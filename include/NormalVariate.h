
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