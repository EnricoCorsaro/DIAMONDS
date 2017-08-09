#include "LinearModel.h"


// LinearModel::LinearModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//

LinearModel::LinearModel(const RefArrayXd covariates)
: Model(covariates)
{
}










// LinearModel::LinearModel()
//
// PURPOSE: 
//      Destructor.
//

LinearModel::~LinearModel()
{

}










// LinearModel::predict()
//
// PURPOSE:
//      Builds the predictions from a linear model of the type f = a*x + b
//      where a and b are the free parameters and x the covariates.
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//      modelParameters:    one-dimensional array where each element
//                          contains the value of a free parameter of the model
//
// OUTPUT:
//      void
//

void LinearModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    // Initialize model parameters

    double slope = modelParameters(0);
    double offset = modelParameters(1);

    
    // Compute predictions

    predictions = slope*covariates + offset;
}

    










