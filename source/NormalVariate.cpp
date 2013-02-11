#include "NormalVariate.h"


// NormalVariate::NormalVariate()
// PURPOSE: 
//      Constructor. Uniform distribution of prior mass (0,1) 
//      and seed number for Marsenne Twister pseudo-random generator are set.
// INPUT:
//      mu: mean of the gaussian distribution
//      sigma: standard deviation of the gaussian distribution
// OUTPUT:

NormalVariate::NormalVariate(double mu, double sigma)
: RandomVariate(1), uniform(0.0, 1.0), engine(time(0)), mu(mu), sigma(sigma)
{

}






// NormalVariate::~NormalVariate()
// PURPOSE: 
//      Destructor
// INPUT:
// OUTPUT:

NormalVariate::~NormalVariate()
{

}






// NormalVariate::drawNestedValues()
// PURPOSE:
//      Draw parameter values from a uniform (flat) proper prior and a one-dimensional Gaussian Likelihood
// INPUT:
//      values: vector to be initialized from prior mass distribution 
//      logLikelihood: vector to be initialized with corresponding logLikelihood values
//      Nvalues: number of objects from nested loop to be initialized
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
//      Draw parameter values from a uniform (flat) proper prior and a one-dimensional Gaussian Likelihood
//      by using the contraint logDensity > logDensityConstraint
//
// INPUT:
//      values: Array to be initialized from prior mass distribution 
//      logDensity: Array to be initialized with corresponding logDensity values
//      logDensityConstraint: constraining value for the logDensity function
//      Nvalues: number of objects from nested loop to be initialized
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
