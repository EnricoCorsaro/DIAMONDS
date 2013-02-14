#include "Prior.h"

// Prior::Prior()
//
// PURPOSE:
//      Constructor.
//
// INPUT:
//
// OUTPUT:

Prior::Prior(int Ndim)
: Ndim(Ndim), 
  uniform(0.0,1.0), 
  engine(time(0))
{
    cerr << "Set parameter space of " << Ndim " dimensions." << endl;
}





// Prior::getNdim()
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:

int Prior::getNdim()
{
    return Ndim;
}





// Prior::setBoundaries()
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:

void Prior::setBoundaries(const RefArrayXd min, const RefArrayXd max)
{
    minimum = min;
    maximum = max;

    return;
}







// Prior::getMinimum()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

double Prior::getMinimum(const int parameterID)
{
    if (parameterID >= minimum.size())
    {
        cerr << "Invalid parameter number. Quitting program.\n";
        exit(1);
    }

    return minimum(parameterID);    
}







// Prior::getMaximum()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

double Prior::getMaximum(const int parameterID)
{
    if (parameterID >= maximum.size())
    {
        cerr << "Invalid parameter number. Quitting program.\n";
        exit(1);
    }

    return maximum(parameterID);    
}








// Prior::uniformPrior()
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:

void Prior::uniformPrior();
{
    uniformFactor = MathExtra::product(1./(maximum - minimum));
    
    return;
}








// Prior::getUniformFactor()
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:

double Prior::getUniformFactor()
{
    return uniformFactor;
}







// Prior::drawFromUniformPrior()
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:
//      void

void Prior::drawFromUniformPrior(RefArrayXXd NestedSampleParameters, const double Nobjects)
{
    ArrayXXd uniformNumbers(Ndim, Nobjects);
    
    for (ptrdiff_t i = 0; i < Ndim; i++)
    {
        for (ptrdiff_t j = 0; j < Nobjects; j++)
        {
            uniformNumbers(i,j) = uniform(engine);
        }
    }

    // TO BE COMPLETED
    return;
}







// Prior::drawFromUniformPriorWithConstraint()
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:

void Prior::drawFromUniformPriorWithConstraint(RefArrayXd NestedSamplerParameter, const double logLikelihoodConstraint)
{
    
    do
    {
        NestedSampleParameter = minimum + uniform(engine)/uniformFactor;
        logLikelihood = MathExtra::logGaussProfile(parameter, mu, sigma);       // TO BE FIXED!!!!!!!
    }
    while (logLikelihood < logLikelihoodConstraint);
    
    return;
}


