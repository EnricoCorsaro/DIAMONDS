#include "ExponentialLikelihood.h"


// ExponentialLikelihood::ExponentialLikelihood()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
// 

ExponentialLikelihood::ExponentialLikelihood(const RefArrayXd observations, Model &model)
: Likelihood(observations, model)
{
}









// ExponentialLikelihood::~ExponentialLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

ExponentialLikelihood::~ExponentialLikelihood()
{
}








// ExponentialLikelihood::logValue()
//
// PURPOSE:
//      Compute the natural logarithm of the exponential likelihood, which
//      s a generalized chi-square distribution with 2 d.o.f.
//
// INPUT:
//      modelParameters: a one-dimensional array containing the actual
//      values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the natural logarithm of the
//      exponential likelihood
//

double ExponentialLikelihood::logValue(RefArrayXd modelParameters)
{
    ArrayXd predictions;
    ArrayXd lambda;
    
    predictions.resize(observations.size());
    model.predict(predictions, modelParameters);
    
    lambda = -1.0*(predictions.log() + observations/predictions);
    
    return lambda.sum();
}

