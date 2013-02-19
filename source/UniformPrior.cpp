
#include "UniformPrior.h"



// UniformPrior:UniformPrior()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
//      minimum: array containing minimum values for setting 
//      lower bounds of the free parameters. 
//      maximum: array containing maximum values for setting 
//      upper bounds of the free parameters.
//      Nobjects: number of objects used for nested sampling.
// 

UniformPrior::UniformPrior(const RefArrayXd min, const RefArrayXd max, const int Nobjects)
: Prior(min.size(), Nobjects), 
  minimum(minimum), 
  maximum(maximum),
  uniform(0.0,1.0),
  engine(time(0))  
{
    assert (min.size() == max.size());

    if (minimum >= maximum)
    {
        cerr << "Invalid boundaries values. Quitting program." << endl;
        exit(1);
    }
    
    uniformFactor = MathExtra::product(1./(maximum - minimum));

    cerr << "Set parameter space of " << Ndimensions " dimensions." << endl;;
} // END UniformPrior::UniformPrior()







// UniformPrior:~UniformPrior()
//
// PURPOSE: 
//      Derived class destructor.
//

UniformPrior::~UniformPrior()
{

} // END UniformPrior::~UniformPrior()







// UniformPrior::getMinimum()
//
// PURPOSE:
//      Get the private data member minimum.
//
// OUTPUT:
//      An array containing the minimum values of the free parameters.
//

ArrayXd UniformPrior::getMinimum()
{
    return minimum;    
} // END UniformPrior::getMinimum()







// UniformPrior::getMaximum()
//
// PURPOSE:
//      Get the private data member maximum.
//
// OUTPUT:
//      An array containing the maximum values of the free parameters.
//

ArrayXd UniformPrior::getMaximum()
{
    return maximum;    
} // END UniformPrior::getMaximum()








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
//
// OUTPUT:
//      void
//

void UniformPrior::draw(RefArrayXXd nestedParameters)
{
    nestedParameters.resize(Ndimensions, Nobjects);
   
    // Uniform sampling over parameters intervals
    for (ptrdiff_t i = 0; i < Ndimensions; i++)
    {
        for (ptrdiff_t j = 0; j < Nobjects; j++)
        {
            nestedParameters(i,j) = uniform(engine)*(maximum(i)-minimum(i)) + minimum(i);
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

void UniformPrior::drawWithConstraint(RefArrayXd nestedParameters, Likelihood &likelihood)
{
    double logLikelihoodConstraint;

    logLikelihoodConstraint = likelihood.logValue(nestedParameters);

    // Uniform sampling to find new parameter with logLikelihood > logLikelihoodConstraint
    do
    {
    
    for (ptrdiff_t i = 0; i < Ndimensions; i++)
        {
            nestedParameters(i) = uniform(engine)*(maximum(i) - minimum(i)) + minimum(i);
        }
    
    logLikelihood = likelihood.logValue(nestedParameters);
    }
    while (logLikelihood < logLikelihoodConstraint);
    
} // END UniformPrior::drawWithConstraint()



