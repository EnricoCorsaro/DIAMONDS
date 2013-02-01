
#include "NormalVariate.h"



// NormalVariate::NormalVariate()
//
// PURPOSE: 
//      Constructor
//
// INPUT:
//      mu: mean of the gaussian distribution
//      sigma: standard deviation of the gaussian distribution
//
// OUTPUT:

NormalVariate::NormalVariate(double mu, double sigma)
: uniform(0.0, 1.0), engine(time(0)), mu(mu), sigma(sigma)
{

}






// NormalVariate::~NormalVariate()
//
// PURPOSE: 
//      Destructor
//
// INPUT:
//
// OUTPUT:

NormalVariate::~NormalVariate()
{

}








// NormalVariate::drawNestedValues()
//
// PURPOSE: 
//
// INPUT:
//      values: 
//      logdensities: 
//      Nvalues:
//
// OUTPUT:

void NormalVariate::drawNestedValues(vector<double> &values, vector<double> &logDensities, int Nvalues)
{
    srand(time(0));
    
    for (int i = 0; i < Nvalues; i++)
    {
        values.at(i) = minimum + uniform(engine) * (maximum - minimum);
    }
    
    MathExtra::logGaussProfile(logDensities, values, mu, sigma);
    return;
}










// NormalVariate::drawNestedValueWithConstraint()
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:

void NormalVariate::drawNestedValueWithConstraint(double &value, double &logDensity, double logDensityConstraint)
{
    srand(time(0));

    do
    {
        value = minimum + uniform(engine) * (maximum - minimum);
        logDensity = MathExtra::logGaussProfile(value, mu, sigma);
    }
    while (logDensity < logDensityConstraint);
    
    return;
}
