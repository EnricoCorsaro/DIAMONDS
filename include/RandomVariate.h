
#ifndef RANDOMVARIATE_H
#define RANDOMVARIATE_H

#include <vector>


using namespace std;

class RandomVariate
{
    public:
    
        RandomVariate();
        ~RandomVariate();
        void setBoundaries(double min, double max);
        virtual void drawNestedValues(vector<double> &values, vector<double> &logDensities, int Nvalues) = 0;
        virtual void drawNestedValueWithConstraint(double &value, double &logDensity, double logDensityConstraint) = 0;
        
    protected:

        double minimum;
        double maximum;
    
    private:

};


#endif