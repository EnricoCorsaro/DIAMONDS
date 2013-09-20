#include "ZeroModel.h"


// ZeroModel::ZeroModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates: one-dimensional array containing the values 
//                  of the independent variable. These will be ignored
//                  as the zero model always returns zero.
//

ZeroModel::ZeroModel(const RefArrayXd covariates)
: Model(covariates)
{

}







// ZeroModel::ZeroModel()
//
// PURPOSE: 
//      Destructor.
//

ZeroModel::~ZeroModel()
{

}







// ZeroModel::predict()
//
// PURPOSE:
//      The predicted values of the ZeroModel are all zero.
//
// INPUT:
//      predictions: One-dimensional array to contain the predictions from the model
//                   These will all be set to zero
//      modelParameters: one-dimensional array containing the value of the free fit parameters
//
// OUTPUT:
//      void
//

void ZeroModel::predict(RefArrayXd predictions, const RefArrayXd modelParameters)
{

    // Check if the predictions array has the right size. If not, resize.

    if (predictions.size() != covariates.size())
    {
        predictions.resize(covariates.size());
    }

    // Fill the predictions array with zeros

    predictions.fill(0.0);

    // That's it.

    return;
}



