
#include "NormalLikelihood.h"


// NormalLikelihood::NormalLikelihood()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
//      covariates: array containing the independent variable values
//      observations: array containing the dependent variable values
//      uncertainties: array containing the uncertainties of the observations
//      model: object specifying the model to be used.
// 

NormalLikelihood::NormalLikelihood(const RefArrayXd covariates, const RefArrayXd observations, const RefArrayXd uncertainties, Model &model)
: Likelihood(covariates, observations, uncertainties, model)
{
    if (covariates.size() != observations.size() || covariates.size() != uncertainties.size())
    {
        cerr << "Array dimensions do not match. Quitting program." << endl;
        exit(EXIT_FAILURE);
    }
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
//      modelParameters: a one-dimensional array containing the actual
//      values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the natural logarithm of the
//      normal likelihood
//

double NormalLikelihood::logValue(RefArrayXd modelParameters)
{
    ArrayXd predictions;
    ArrayXd lambda;
    ArrayXd lambda0;
    
    predictions.resize(covariates.size());
    model.predict(predictions, modelParameters);
    
    lambda0 = -0.5 * observations.size() * log(2.0*MathExtra::PI) -1.0 * uncertainties.log(); 
    lambda = lambda0 - 0.5 * ((observations - predictions)*(observations - predictions)) / (uncertainties*uncertainties);
    
    return lambda.sum();
} // END NormalLikelihood::logValue()



