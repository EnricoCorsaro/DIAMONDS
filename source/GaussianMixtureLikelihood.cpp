#include "GaussianMixtureLikelihood.h"


// GaussianMixtureLikelihood::GaussianMixtureLikelihood()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
//      observations: array containing the dependent variable values
//      uncertainties: array containing the uncertainties of the observations
//      model: object specifying the model to be used.
// 

GaussianMixtureLikelihood::GaussianMixtureLikelihood(const RefArrayXd observations, const RefArrayXd uncertainties, Model &model)
: Likelihood(observations, model),
  uncertainties(uncertainties)
{
    assert(observations.size() == uncertainties.size());
}









// GaussianMixtureLikelihood::~GaussianMixtureLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

GaussianMixtureLikelihood::~GaussianMixtureLikelihood()
{

}









// GaussianMixtureLikelihood::getUncertainties();
//
// PURPOSE:
//      Get protected data member uncertainties.
//
// OUTPUT:
//      uncertainties: one-dimensional array containing the
//      uncertainties on the dependent variable values.
//

ArrayXd GaussianMixtureLikelihood::getUncertainties()
{
    return uncertainties;
}








// GaussianMixtureLikelihood::logValue()
//
// PURPOSE:
//      Compute the natural logarithm of the normal likelihood for 
//      a given set of observations, uncertainties and predictions.
//
// INPUT:
//      modelParameters: a one-dimensional array containing the actual
//      values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the natural logarithm of the
//      normal likelihood
//

double GaussianMixtureLikelihood::logValue(RefArrayXd modelParameters)
{
    ArrayXd predictions;
    ArrayXd lambda;
    ArrayXd argument1;
    ArrayXd argument2;
    double weight = modelParameters(4);

    int Npoints = observations.size();
    predictions.resize(Npoints*2);
    predictions.setZero();
    lambda.resize(Npoints);
    model.predict(predictions, modelParameters);
   
    argument1 = -0.5 * ((observations - predictions.segment(0,Npoints))*(observations - predictions.segment(0,Npoints))) / (uncertainties*uncertainties);
    argument2 = -0.5 * ((observations - predictions.segment(Npoints-1,Npoints))*(observations - predictions.segment(Npoints-1,Npoints))) / (uncertainties*uncertainties);
    
    for (int i = 0; i < Npoints; i++)
    {
        if (argument1(i) >= argument2(i))
        {
            lambda(i) = -0.5 * log(2.0*Functions::PI) -1.0 * log(uncertainties(i)) + log(weight) + argument1(i) + log(1.0 + ((1.0 - weight)/weight) * exp(argument2(i) - argument1(i)));
        }
        else
        {
            lambda(i) = -0.5 * log(2.0*Functions::PI) -1.0 * log(uncertainties(i)) + log(1.0 - weight) + argument2(i) + log(1.0 + (weight/(1.0 - weight)) * exp(argument1(i) - argument2(i)));
        }
    }

    return lambda.sum();
}



