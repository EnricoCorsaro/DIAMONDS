#include "NormalLikelihood.h"


// NormalLikelihood::NormalLikelihood()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
//      observations: array containing the dependent variable values
//      uncertainties: array containing the uncertainties of the observations
//      model: object specifying the model to be used.
// 

NormalLikelihood::NormalLikelihood(const RefArrayXd observations, const RefArrayXd uncertainties, Model &model)
: Likelihood(observations, uncertainties, model)
{
    assert(observations.size() || uncertainties.size());
} // END NormalLikelihood::NormalLikelihood()









// NormalLikelihood::~NormalLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

NormalLikelihood::~NormalLikelihood()
{

} // END NormalLikelihood::~NormalLikelihood()









// NormalLikelihood::getUncertainties();
//
// PURPOSE:
//      Get protected data member uncertainties.
//
// OUTPUT:
//      uncertainties: one-dimensional array containing the
//      uncertainties on the dependent variable values.
//

ArrayXd NormalLikelihood::getUncertainties()
{
    return uncertainties;
} // END NormalLikelihood::getUncertainties()








// NormalLikelihood::logValue()
//
// PURPOSE:
//      Compute the natural logarithm of the normal likelihood for 
//      a given set of observations, uncertainties and predictions.
//
// INPUT:
//      nestedSampleOfParameters: a one-dimensional array containing the actual
//      values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the natural logarithm of the
//      normal likelihood
//

double NormalLikelihood::logValue(RefArrayXd nestedSampleOfParameters)
{
    ArrayXd predictions;
    ArrayXd lambda;
    ArrayXd lambda0;
    
    predictions.resize(observations.size());
    model.predict(predictions, nestedSampleOfParameters);
    
    lambda0 = -0.5 * observations.size() * log(2.0*Functions::PI) -1.0 * uncertainties.log(); 
    lambda = lambda0 - 0.5 * ((observations - predictions)*(observations - predictions)) / (uncertainties*uncertainties);
    
    return lambda.sum();
} // END NormalLikelihood::logValue()



