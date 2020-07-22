#include "MultiLinearNormalLikelihood.h"


// MultiLinearNormalLikelihood::MultiLinearNormalLikelihood()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
//      observations: array containing the dependent variable values
//      uncertainties: array containing the uncertainties of the observations
//      model: object specifying the model to be used.
// 

MultiLinearNormalLikelihood::MultiLinearNormalLikelihood(const RefArrayXd observations, const RefArrayXd covariatesUncertainties, 
                                                         const RefArrayXd observationsUncertainty, MultiLinearModel &model)
: Likelihood(observations, model),
  covariatesUncertainties(covariatesUncertainties),
  observationsUncertainty(observationsUncertainty)
{
    Nobservables = model.getNobservables();
    Npoints = model.getNpoints();
    assert(observations.size() == observationsUncertainty.size());
}









// MultiLinearNormalLikelihood::~MultiLinearNormalLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

MultiLinearNormalLikelihood::~MultiLinearNormalLikelihood()
{

}









// MultiLinearNormalLikelihood::getCovariatesUncertainties();
//
// PURPOSE:
//      Get protected data member covariatesUncertainties.
//
// OUTPUT:
//      covariatesUncertainties: one-dimensional array containing the
//      uncertainties on the independent variables.
//

ArrayXd MultiLinearNormalLikelihood::getCovariatesUncertainties()
{
    return covariatesUncertainties;
}








// MultiLinearNormalLikelihood::getObservationsUncertainty();
//
// PURPOSE:
//      Get protected data member observationsUncertainty.
//
// OUTPUT:
//      observationsUncertainty: one-dimensional array containing the
//      uncertainties on the dependent variable.
//

ArrayXd MultiLinearNormalLikelihood::getObservationsUncertainty()
{
    return observationsUncertainty;
}











// MultiLinearNormalLikelihood::logValue()
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
//      For this multilinear version of the normal likelihood, it is required
//      that the last element of the modelParameters is the free parameter of the offset.
//
// OUTPUT:
//      a double number containing the natural logarithm of the
//      normal likelihood
//

double MultiLinearNormalLikelihood::logValue(RefArrayXd modelParameters)
{
    ArrayXd predictions;
    ArrayXd lambda;
    ArrayXd lambda0;
    ArrayXd totalVariance(Npoints);

    predictions.resize(observations.size());
    predictions.setZero();
    model.predict(predictions, modelParameters);

    totalVariance = observationsUncertainty.square();

    for (int observable = 0; observable < Nobservables; ++observable)
    {
        totalVariance += covariatesUncertainties.segment(observable*Npoints, Npoints).square()
                        *modelParameters(observable)*modelParameters(observable);
    }
    
    lambda0 = -0.5 * log(2.0*Functions::PI) -0.5 * totalVariance.log(); 
    lambda = lambda0 - 0.5 * ((observations - predictions)*(observations - predictions)) / totalVariance;
    
    return lambda.sum();
}



