
#include "NormalLikelihood.h"


// NormalLikelihood::NormalLikelihood()
//
// PURPOSE: constructor
//
// INPUT:
//      covariates:
//      observations:
//      uncertainties:
//      model:
// 

NormalLikelihood::NormalLikelihood(RefArrayXd covariates, RefArrayXd observations, RefArrayXd uncertainties, Model &model)
: Likelihood(covariates, observations, uncertainties, model)
{

}





// NormalLikelihood::~NormalLikelihood()
//
// PURPOSE: destructor
//

NormalLikelihood::~NormalLikelihood()
{

}






// NormalLikelihood::~NormalLikelihood()
//
// PURPOSE: destructor
//

double NormalLikelihood::logDensity(RefArrayXd modelParameters)
{
    ArrayXd predictions;
    model.predict(predictions, modelParameters);
    
    double term = -0.5 * observations.size() * cmath.log(2.0*MathExtra::PI); 
    delta = -uncertainties - 0.5 * ((observations - predictions)*(observations - predictions)) / (uncertainties*uncertainties);
    return term + delta.sum();
}