
#ifndef RANDOMVARIATE_H
#define RANDOMVARIATE_H

#include <cassert>
#include <Eigen/Core>

using namespace std;

using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class RandomVariate
{
    public:
    
        RandomVariate(int Ndim);
        ~RandomVariate();
        int getNdim();
        void setBoundaries(const RefArrayXd min, const RefArrayXd max);
        void setBoundaries(double min, double max);
        
        virtual void drawNestedValues(RefArrayXXd values, RefArrayXd logDensities, int Nvalues) = 0;
        virtual void drawNestedValueWithConstraint(RefArrayXd value, double &logDensity, double logDensityConstraint) = 0;
        
    protected:

        int Ndim;
        ArrayXd minimum;
        ArrayXd maximum;
    
    private:

};



#endif