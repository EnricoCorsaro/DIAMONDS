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

double ExponentialLikelihood::logValue(RefArrayXd const modelParameters)
{
    ArrayXd predictions(observations.size());
    predictions.setZero();
    model.predict(predictions, modelParameters);

    return -1.0*(predictions.log() + observations/predictions).sum();
}

