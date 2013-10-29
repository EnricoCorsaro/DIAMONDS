#include "DeltaPrior.h"



// DeltaPrior::DeltaPrior()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
//      constantParameters:  an eigen array containing the fixed parameters to be
//                          considered as constants within the computation.
//

DeltaPrior::DeltaPrior(const RefArrayXd constantParameters)
: Prior(constantParameters.size()),
  constantParameters(constantParameters)
{
}









// DeltaPrior:~DeltaPrior()
//
// PURPOSE: 
//      Derived class destructor.
//

DeltaPrior::~DeltaPrior()
{

}










// DeltaPrior::getConstant()
//
// PURPOSE:
//      Get the private data member constantParameters.
//
// OUTPUT:
//      An eigen array containing the value of the constant parameters.
//

ArrayXd DeltaPrior::getConstants()
{
    return constantParameters;    
}










// DeltaPrior::logDensity()
//
// PURPOSE:
//      Compute the logarithm of the probability density distribution evaluated in 'x'.
//
// INPUT: 
//      x:                      Point in which the log(pdf) should be evaluated.
//      includeConstantTerm:    If true : compute the exact log(density), 
//                              If false: ignore the constant terms (with factors of pi, 2, etc.)
//
// OUTPUT:
//      Natural logarithm of the probability density evaluation in x.
//

double DeltaPrior::logDensity(RefArrayXd x, const bool includeConstantTerm)
{
    double logDens;

    if ((x != constantParameters).any())
    {
        // The point x is not the one allowed by the delta prior. In this case
        // the density is zero, and thus the log(density) is -infinity. The value
        // of the latter is defined in the parent class.

        logDens = minusInfinity;
        return logDens;
    }
    else
    {
        // The point is exactly at the only possible set of coordinates.

        logDens = 0.0;
    }

    return logDens;
}















// DeltaPrior::draw()
//
// PURPOSE:
//      Draw the sample of parameters values from the delta prior
//      distribution. The parameters are in number Ndimensions
//      and contain Npoints values each.
//
// INPUT:
//      drawnSample:    two-dimensional Eigen Array to contain 
//                      the resulting parameters values.
//
// OUTPUT:
//      void
//

void DeltaPrior::draw(RefArrayXXd drawnSample)
{
    int Npoints = drawnSample.cols();
 

    // Delta sampling over all free parameters and points

    for (int i = 0; i < Ndimensions; i++)
    {
        for (int j = 0; j < Npoints; j++)
        {
            drawnSample(i,j) = constantParameters(i);
        }
    }
}
