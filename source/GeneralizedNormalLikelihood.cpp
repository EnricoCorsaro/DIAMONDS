#include "GeneralizedNormalLikelihood.h"


// GeneralizedNormalLikelihood::GeneralizedNormalLikelihood()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
//      observations: array containing the dependent variable values
//      observationsUncertainties: array containing the uncertainties of the observations
//      model: object specifying the model to be used.
// 

GeneralizedNormalLikelihood::GeneralizedNormalLikelihood(const RefArrayXd observations, const RefArrayXd observationsUncertainty, 
                                                        Model &model)
: Likelihood(observations, model),
  observationsUncertainty(observationsUncertainty)
{
    assert(observations.size() == observationsUncertainty.size());
}









// GeneralizedNormalLikelihood::~GeneralizedNormalLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

GeneralizedNormalLikelihood::~GeneralizedNormalLikelihood()
{

}











// GeneralizedNormalLikelihood::getObservationsUncertainty();
//
// PURPOSE:
//      Get protected data member observationsUncertainty.
//
// OUTPUT:
//      observationsUncertainty: one-dimensional array containing the
//      uncertainties on the dependent variable.
//

ArrayXd GeneralizedNormalLikelihood::getObservationsUncertainty()
{
    return observationsUncertainty;
}











// GeneralizedNormalLikelihood::logValue()
//
// PURPOSE:
//      Compute the natural logarithm of the normal likelihood for 
//      a given set of observations, uncertainties and predictions of a
//      multilinear model. In this case, also uncertainties depend upon
//      the free parameters of the model.
//
// INPUT:
//      modelParameters: a one-dimensional array containing the actual
//      values of the free parameters that describe the model.
//      For this version of the normal likelihood, it is required
//      that the last element of the modelParameters is the free parameter of the offset and
//      that the remainder of the elements are divided in two groups, the first one related to
//      the multi-linear part of the model, and the second one related to the polynomial part of the model.
//
// OUTPUT:
//      a double number containing the natural logarithm of the
//      normal likelihood
//

double GeneralizedNormalLikelihood::logValue(RefArrayXd modelParameters)
{
    ArrayXd predictions;
    ArrayXd totalVariance;
    ArrayXd lambda;
    ArrayXd lambda0;

    predictions.resize(observations.size());
    predictions.setZero();
    totalVariance.resize(observations.size());
    totalVariance.setZero();

    model.predict(predictions, modelParameters);
    model.computeVariance(totalVariance, modelParameters);

    totalVariance += observationsUncertainty.square();

    lambda0 = -0.5 * log(2.0*Functions::PI) -0.5 * totalVariance.log(); 
    lambda = lambda0 - 0.5 * ((observations - predictions)*(observations - predictions)) / totalVariance;
    
    return lambda.sum();
}



