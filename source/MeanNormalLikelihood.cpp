
#include "MeanNormalLikelihood.h"


// MeanNormalLikelihood::MeanNormalLikelihood()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
//      observations: array containing the dependent variable values
//      uncertainties: array containing the uncertainties of the observations
//      model: object specifying the model to be used.
// 

MeanNormalLikelihood::MeanNormalLikelihood(const RefArrayXd observations, const RefArrayXd uncertainties, Model &model)
: Likelihood(observations, model),
  uncertainties(uncertainties)
{
    double normalizeFactor;
    assert(observations.size() == uncertainties.size());
    normalizeFactor = sqrt(observations.size()/uncertainties.pow(-2).sum());
    normalizedUncertainties = uncertainties/normalizeFactor; 
    weights = normalizedUncertainties.pow(-2);

} // END MeanNormalLikelihood::MeanNormalLikelihood()








// MeanNormalLikelihood::~MeanNormalLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

MeanNormalLikelihood::~MeanNormalLikelihood()
{

} // END MeanNormalLikelihood::~MeanNormalLikelihood()










// MeanNormalLikelihood::getNormalizedUncertainties();
//
// PURPOSE:
//      Get private data member normalizedUncertainties.
//
// OUTPUT:
//      uncertainties: one-dimensional array containing the
//      uncertainties on the dependent variable values.
//

ArrayXd MeanNormalLikelihood::getNormalizedUncertainties()
{
    return normalizedUncertainties;
} // END MeanNormalLikelihood::getNormalizedUncertainties()








// MeanNormalLikelihood::logValue()
//
// PURPOSE:
//      Compute the natural logarithm of the mean normal likelihood for 
//      a given set of observations, uncertainties and predictions.
//      The mean normal likelihood is a normal likelihood whose uncertainties
//      have been integrated away by means of a Jeffrey's prior.
//      For more details cf. Froehlich H.-E. et al. 2009, A&A, 506, 263.
//
// INPUT:
//      modelParameters: a one-dimensional array containing the actual
//      values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the natural logarithm of the
//      mean likelihood
//

double MeanNormalLikelihood::logValue(RefArrayXd modelParameters)
{
    unsigned long n = observations.size();
    double lambda0;
    double lambda;
    ArrayXd argument;
    ArrayXd predictions;

    predictions.resize(n);
    predictions.setZero();
    model.predict(predictions, modelParameters);
    argument = (observations - predictions);
    argument = argument.square()*weights;

    lambda0 = lgammal(n/2.) - log(2) - (n/2.)*log(Functions::PI) + 0.5*weights.log().sum();
    lambda = lambda0 - (n/2.)*log(argument.sum());

    return lambda;
}



