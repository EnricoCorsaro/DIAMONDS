
#include "Likelihood.h"



// Likelihood::Likelihood()
//
// PURPOSE: 
//      Abstract base class constructor.
//
// INPUT:
//      covariates: array containing the independent variable values
//      observations: array containing the dependent variable values
//      uncertainties: array containing the uncertainties of the observations
//      model: object specifying the model to be used.
// 

Likelihood::Likelihood(const RefArrayXd covariates, const RefArrayXd observations, const RefArrayXd uncertainties, Model &model)
: covariates(covariates), 
  observations(observations), 
  uncertainties(uncertainties), 
  model(model)
{

} // END Likelihood::Likelihood()









// Likelihood::~Likelihood()
//
// PURPOSE: 
//      Abstract base class destructor.
//

Likelihood::~Likelihood()
{

} // END Likelihood::Likelihood()










// Likelihood::getCovariates();
//
// PURPOSE:
//      Get protected data member covariates.
//
// OUTPUT:
//      covariates: one-dimensional array containing the
//      independent variable values.
//

ArrayXd Likelihood::getCovariates()
{
    return covariates;
} // END Likelihood::getCovariates()









// Likelihood::getObservations();
//
// PURPOSE:
//      Get protected data member observations.
//
// OUTPUT:
//      observations: one-dimensional array containing the
//      dependent variable values.
//

ArrayXd Likelihood::getObservations()
{
    return observations;
} // END Likelihood::getObservations()



