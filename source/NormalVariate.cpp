
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
: RandomVariate(1), uniform(0.0, 1.0), engine(time(0)), mu(mu), sigma(sigma)
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


void NormalVariate::drawNestedValues(RefArrayXXd values, RefArrayXd logDensities, int Nvalues)
{
    ArrayXd uniformNumbers(Nvalues);
    for (int i = 0; i < Nvalues; i++) 
    {
        uniformNumbers(i) = uniform(engine);
    }

    values.row(0) = minimum(0) + uniformNumbers * (maximum(0) - minimum(0));
    
    MathExtra::logGaussProfile(logDensities, values.row(0), mu, sigma, 1.0);
    return;
}










// NormalVariate::drawNestedValueWithConstraint()
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:

void NormalVariate::drawNestedValueWithConstraint(RefArrayXd value, double &logDensity, double logDensityConstraint)
{
    do
    {
        value(0) = minimum(0) + uniform(engine) * (maximum(0) - minimum(0));
        logDensity = MathExtra::logGaussProfile(value(0), mu, sigma);
    }
    while (logDensity < logDensityConstraint);
    
    return;
}
