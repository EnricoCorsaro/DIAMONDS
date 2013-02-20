
#include "UniformPrior.h"



// UniformPrior:UniformPrior()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
//      boundaries: array containing minimum and maximum values for setting 
//      lower and upper bounds of the free parameters. 
//      Nobjects: number of objects used for nested sampling.
// 
// NOTE:
//      boundaries array is a (Ndimensions * 2) format matrix, where the first
//      column contains the minimum values and the second column the corresponding
//      maximum values for each free parameter.
//

UniformPrior::UniformPrior(const RefArrayXXd boundaries, const int Nobjects)
: Prior(boundaries.size()/2, Nobjects), 
  uniform(0.0,1.0),
  engine(time(0)),  
  boundaries(boundaries)  
{
    uniformFactor = (1./(boundaries.col(1) - boundaries.col(0))).prod();

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











// UniformPrior::getBoundaries()
//
// PURPOSE:
//      Get the private data member boundaries.
//
// OUTPUT:
//      An array containing the minimum and maximum values of the free parameters.
//

ArrayXXd UniformPrior::getBoundaries()
{
    return boundaries;    
} // END UniformPrior::getBoundaries()









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
    // Uniform sampling over parameters intervals
    for (ptrdiff_t i = 0; i < Ndimensions; i++)
    {
        for (ptrdiff_t j = 0; j < Nobjects; j++)
        {
            nestedParameters(i,j) = uniform(engine)*(boundaries(i,1)-boundaries(i,0)) + boundaries(i,0);
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
        //for (ptrdiff_t i = 0; i < Ndimensions; i++)
        //    {
        //        nestedParameters(i) = uniform(engine)*(boundaries(i,1) - boundaries(i,0)) + boundaries(i,0);
        //    }
    
        nestedParameters = uniform(engine)*(boundaries.col(1) - boundaries.col(0)) + boundaries.col(0);
        logLikelihood = likelihood.logValue(nestedParameters);
    }
    while (logLikelihood < logLikelihoodConstraint);
    
} // END UniformPrior::drawWithConstraint()


