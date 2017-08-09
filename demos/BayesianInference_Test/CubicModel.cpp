#include "CubicModel.h"


// CubicModel::CubicModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//

CubicModel::CubicModel(const RefArrayXd covariates)
: Model(covariates)
{
}










// CubicModel::CubicModel()
//
// PURPOSE: 
//      Destructor.
//

CubicModel::~CubicModel()
{

}










// CubicModel::predict()
//
// PURPOSE:
//      Builds the predictions from a quadratic model of the type f = a*x + b*x^2 + c*x^3 + d
//      where a, b, c, and d are the free parameters and x the covariates.
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

void CubicModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    // Initialize model parameters

    double linearSlope = modelParameters(0);
    double quadraticSlope = modelParameters(1);
    double cubicSlope = modelParameters(2);
    double offset = modelParameters(3);

    
    // Compute predictions

    predictions = linearSlope*covariates + quadraticSlope*covariates*covariates + cubicSlope*covariates*covariates*covariates + offset;
}

    










