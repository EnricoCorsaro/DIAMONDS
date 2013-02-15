
#include "Likelihood.h"



// Likelihood::Likelihood()
//
// PURPOSE: constructor
//
// INPUT:
//      covariates:
//      observations:
//      uncertainties:
//      model:
// 

Likelihood::Likelihood(RefArrayXd covariates, RefArrayXd observations, RefArrayXd uncertainties, Model &model);
: covariates(covariates), observations(observations), uncertainties(uncertainties), model(model)
{

}





// Likelihood::~Likelihood()
//
// PURPOSE: destructor
//

Likelihood::~Likelihood()
{

}

