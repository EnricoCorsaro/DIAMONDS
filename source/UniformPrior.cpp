
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
  engine(time(0)),  
  minima(minima),
  maxima(maxima)
{
    assert (minima.size() == maxima.size());
    uniformFactor = (1./(maxima - minima)).prod();
    cerr << "Set parameter space of " << Ndimensions << " dimensions." << endl;
} // END UniformPrior::UniformPrior()









// UniformPrior:~UniformPrior()
//
// PURPOSE: 
//      Derived class destructor.
//

UniformPrior::~UniformPrior()
{

} // END UniformPrior::~UniformPrior()











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
} // END UniformPrior::getMinima()











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
} // END UniformPrior::getMaxima()










// UniformPrior::getUniformFactor()
//
// PURPOSE: 
//      Get the private data member uniformFactor.
//
// OUTPUT:
//      An integer containing the normalization factor of the uniform prior.

double UniformPrior::getUniformFactor()
{
    return uniformFactor;
} // END UniformPrior::getUniformPrior()








// UniformPrior::draw()
//
// PURPOSE:
//      Draw a sample of parameters values from a uniform prior
//      distributions. The parameters are in number Ndimensions
//      and contain Nobjects values each.
//
// INPUT:
//      nestedParameters: two-dimensional array to contain 
//      the resulting parameters values.
//      Nobjects: integer containing the number of objects used 
//      in the nested sampling process.
//
// OUTPUT:
//      void
//

void UniformPrior::draw(RefArrayXXd nestedParameters, const int Nobjects)
{
    // Uniform sampling over parameters intervals
    
    for (int i = 0; i < Ndimensions; i++)
    {
        for (int j = 0; j < Nobjects; j++)
        {
            nestedParameters(i,j) = uniform(engine)*(maxima(i)-minima(i)) + minima(i);
        }
    }

} // END UniformPrior::draw()







// UniformPrior::drawWithConstraint()
//
// PURPOSE: 
//      Replace an old set of parameters values with a new one
//      having higher likelihood value.
//
// INPUT:
//      nestedParameters: one-dimensional array containing the set of 
//      parameters values to be updated.
//      likelihood: an object to compute the corresponding likelihood value.
//
// OUTPUT:
//      void
//
// NOTE:
//      nestedParameters refers to the worst object identified in the nested
//      sampling loop. Thus, the array contains Ndimensions elements.
//

void UniformPrior::drawWithConstraint(RefArrayXd nestedParameters, Likelihood &likelihood)
{
    double logLikelihood;
    double logLikelihoodConstraint = likelihood.logValue(nestedParameters);
    
    // Uniform sampling to find new parameter with logLikelihood > logLikelihoodConstraint
    
    do
    {
        for (int i = 0; i < Ndimensions; i++)
            {
                nestedParameters(i) = uniform(engine)*(maxima(i) - minima(i)) + minima(i);
            }
    
        // nestedParameters = uniform(engine)*(maxima - minima) + minima;
        
        logLikelihood = likelihood.logValue(nestedParameters);
    }
    while (logLikelihood < logLikelihoodConstraint);
    
} // END UniformPrior::drawWithConstraint()


