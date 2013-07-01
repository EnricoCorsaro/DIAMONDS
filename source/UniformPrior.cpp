#include "UniformPrior.h"



// UniformPrior:UniformPrior()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
//      minima: array containing minimum values for setting 
//      lower bounds of the free parameters. 
//      maxima: array containing maximum values for setting
//      upper bounds of the free parameters.
//

UniformPrior::UniformPrior(const RefArrayXd minima, const RefArrayXd maxima)
: Prior(minima.size()),
  uniform(0.0,1.0),
  minima(minima),
  maxima(maxima)
{
    assert (minima.size() == maxima.size());
}









// UniformPrior:~UniformPrior()
//
// PURPOSE: 
//      Derived class destructor.
//

UniformPrior::~UniformPrior()
{

}









// UniformPrior::getMinima()
//
// PURPOSE:
//      Get the private data member minima.
//
// OUTPUT:
//      An array containing the minimum values of the free parameters.
//

ArrayXd UniformPrior::getMinima()
{
    return minima;    
}











// UniformPrior::getMaxima()
//
// PURPOSE:
//      Get the private data member maxima.
//
// OUTPUT:
//      An array containing the maximum values of the free parameters.
//

ArrayXd UniformPrior::getMaxima()
{
    return maxima;    
}










// UniformPrior::logDensity()
//
// PURPOSE:
//      Compute the logarithm of the probability density distribution evaluated in 'x'.
//
// INPUT: 
//      x:                   Point in which the log(pdf) should be evaluated.
//      includeConstantTerm: If true : compute the exact log(density), 
//                           If false: ignore the constant terms (with factors of pi, 2, etc.)
//
// OUTPUT:
//      Natural logarithm of the probability density evaluation in x.
//

double UniformPrior::logDensity(RefArrayXd x, const bool includeConstantTerm)
{
    double logDens;

    if ((x < minima).any() | (x > maxima).any())
    {
        // The point x falls out of the distribution's boundaries. In this case
        // the density is zero, and thus the log(density) is infinity, which we
        // approximate with the largest possible value.

        logDens = numeric_limits<double>::max();
        return logDens;
    }
    else
    {
        // The points falls inside the distribution's boundaries.

        logDens = -1.0;
    }

    if (includeConstantTerm)
    {
        logDens *= (maxima - minima).log().sum(); 
    }

    return logDens;
}















// UniformPrior::draw()
//
// PURPOSE:
//      Draw a sample of parameters values from a uniform prior
//      distributions. The parameters are in number Ndimensions
//      and contain Nobjects values each.
//
// INPUT:
//      sample: two-dimensional Eigen Array to contain 
//              the resulting parameters values.
//
// OUTPUT:
//      void
//

void UniformPrior::draw(RefArrayXXd sample)
{
    // Uniform sampling over parameters intervals
    
    for (int i = 0; i < Ndimensions; i++)
    {
        for (int j = 0; j < sample.cols(); j++)
        {
            sample(i,j) = uniform(engine)*(maxima(i)-minima(i)) + minima(i);
        }
    }

}







// UniformPrior::drawWithConstraint()
//
// PURPOSE: 
//      Replace an old set of parameters values with a new one
//      having higher likelihood value.
//
// INPUT:
//      parameters: one-dimensional Eigen Array containing the set of 
//      parameters values to be updated.
//      likelihood: an object to compute the corresponding likelihood value.
//
// OUTPUT:
//      void
//
// NOTE:
//      parameters refers to the worst object identified in the nested
//      sampling loop. Thus, the array contains Ndimensions elements.
//

void UniformPrior::drawWithConstraint(RefArrayXd parameters, Likelihood &likelihood)
{
    double logLikelihood;
    double logLikelihoodConstraint = likelihood.logValue(parameters);


    // Uniform sampling to find new parameter with logLikelihood > logLikelihoodConstraint
    
    do
    {
        for (int i = 0; i < Ndimensions; i++)
            {
                parameters(i) = uniform(engine)*(maxima(i) - minima(i)) + minima(i);
            }
    
        logLikelihood = likelihood.logValue(parameters);
    }
    while (logLikelihood <= logLikelihoodConstraint);
    

}



