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
: uniform(0.0, 1.0), engine(time(0)), mu(mu), sigma(sigma)
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

void NormalVariate::drawNestedValues(vector<double> &values, vector<double> &logLikelihood, int Nvalues)
{
    for (int i = 0; i < Nvalues; i++)
    {
        values.at(i) = minimum + uniform(engine) * (maximum - minimum);     // Converting prior mass values to parameter values
    }
    
    MathExtra::logGaussProfile(logLikelihood, values, mu, sigma);           // Evaluate logLikelihood values of the parameters
    return;
}





// NormalVariate::drawNestedValueWithConstraint()
// PURPOSE: 
//      Draw parameter values from a uniform (flat) proper prior and a one-dimensional Gaussian Likelihood
//      by using the contraint logLikelihood > logLikelihood*
// INPUT:
//      values: vector to be initialized from prior mass distribution 
//      logLikelihood: vector to be initialized with corresponding logLikelihood values
//      logLikelihoodConstraint: constraining value for the logLikelihood function
//      Nvalues: number of objects from nested loop to be initialized
// OUTPUT:

void NormalVariate::drawNestedValueWithConstraint(double &value, double &logLikelihood, double logLikelihoodConstraint)
{
    do
    {
        value = minimum + uniform(engine) * (maximum - minimum);
        logLikelihood = MathExtra::logGaussProfile(value, mu, sigma);
    }
    while (logLikelihood < logLikelihoodConstraint);
    
    return;
}
