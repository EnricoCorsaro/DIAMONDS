#include "QuadraticModel.h"


// QuadraticModel::QuadraticModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//

QuadraticModel::QuadraticModel(const RefArrayXd covariates)
: Model(covariates)
{
}










// QuadraticModel::QuadraticModel()
//
// PURPOSE: 
//      Destructor.
//

QuadraticModel::~QuadraticModel()
{

}










// QuadraticModel::predict()
//
// PURPOSE:
//      Builds the predictions from a quadratic model of the type f = a*x + b*x^2 + c
//      where a, b, and c are the free parameters and x the covariates.
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

void QuadraticModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    // Initialize model parameters

    double linearSlope = modelParameters(0);
    double quadraticSlope = modelParameters(1);
    double offset = modelParameters(2);

    
    // Compute predictions

    predictions = linearSlope*covariates + quadraticSlope*covariates*covariates + offset;
}

    










