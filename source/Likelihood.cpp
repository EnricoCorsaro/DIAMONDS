
#include "Likelihood.h"



// Likelihood::Likelihood()
//
// PURPOSE: 
//      Abstract base class constructor.
//
// INPUT:
//      observations: array containing the dependent variable values
//      model: object specifying the model to be used.
// 

Likelihood::Likelihood(const RefArrayXd observations, Model &model)
: observations(observations), 
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



