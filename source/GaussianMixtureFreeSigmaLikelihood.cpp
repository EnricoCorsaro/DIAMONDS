#include "GaussianMixtureFreeSigmaLikelihood.h"


// GaussianMixtureFreeSigmaLikelihood::GaussianMixtureFreeSigmaLikelihood()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
//      observations: array containing the dependent variable values
//      model: object specifying the model to be used.
// 

GaussianMixtureFreeSigmaLikelihood::GaussianMixtureFreeSigmaLikelihood(const RefArrayXd observations, Model &model)
: Likelihood(observations, model)
{
}









// GaussianMixtureFreeSigmaLikelihood::~GaussianMixtureFreeSigmaLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

GaussianMixtureFreeSigmaLikelihood::~GaussianMixtureFreeSigmaLikelihood()
{

}











// GaussianMixtureFreeSigmaLikelihood::logValue()
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

double GaussianMixtureFreeSigmaLikelihood::logValue(RefArrayXd modelParameters)
{
    ArrayXd predictions;
    ArrayXd lambda;
    ArrayXd argument1;
    ArrayXd argument2;
    double weight = modelParameters(4);
    double logSigma = modelParameters(5);

    int Npoints = observations.size();
    
    predictions.resize(Npoints*2);
    predictions.setZero();
    lambda.resize(Npoints);
    
    model.predict(predictions, modelParameters);
    ArrayXd logbit;
 
    argument1 = -0.5 * ((observations - predictions.segment(0,Npoints))*(observations - predictions.segment(0,Npoints))) / (exp(2*logSigma));
    argument2 = -0.5 * ((observations - predictions.segment(Npoints-1,Npoints))*(observations - predictions.segment(Npoints-1,Npoints))) / (exp(2*logSigma));
   
    for (int i = 0; i < Npoints; i++)
    {
        if (argument1(i) >= argument2(i))
        {
            lambda(i) = -0.5 * log(2.0*Functions::PI) -1.0 * logSigma + log(weight) + argument1(i) + log(1.0 + ((1.0 - weight)/weight) * exp(argument2(i) - argument1(i)));
        }
        else
        {
            lambda(i) = -0.5 * log(2.0*Functions::PI) -1.0 * logSigma + log(1.0 - weight) + argument2(i) + log(1.0 + (weight/(1.0 - weight)) * exp(argument1(i) - argument2(i)));
        }
    } 
   
    return lambda.sum();
}



