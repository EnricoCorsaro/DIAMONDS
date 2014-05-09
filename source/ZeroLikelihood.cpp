#include "ZeroLikelihood.h"


// ZeroLikelihood::ZeroLikelihood()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
//      observations: array containing the dependent variable values
//      model: object specifying the model to be used.
// 

ZeroLikelihood::ZeroLikelihood(const RefArrayXd observations, Model &model)
: Likelihood(observations, model)
{
}









// ZeroLikelihood::~ZeroLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

ZeroLikelihood::~ZeroLikelihood()
{

}









// ZeroLikelihood::logValue()
//
// PURPOSE:
//      Returns a zero value for the likelihood, which implies
//      a -infinity for its natural logarithm.
//
// INPUT:
//      modelParameters: a one-dimensional array containing the actual
//      values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the natural logarithm of zero.
//

double ZeroLikelihood::logValue(RefArrayXd modelParameters)
{
    double zeroValue = numeric_limits<double>::lowest();
    return zeroValue;
}



