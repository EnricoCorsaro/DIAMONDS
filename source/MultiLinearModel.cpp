#include "MultiLinearModel.h"


// MultiLinearModel::MultiLinearModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      multiCovariates:        one-dimensional array containing the values
//                              of the independent variables concatenated in the same
//                              order as the columns from left to right in the input data file.
//      Nobservables:           the number of observables in the multilinear model
//

MultiLinearModel::MultiLinearModel(const RefArrayXd covariates, int Nobservables)
: Model(covariates),
  Nobservables(Nobservables)
{
    Npoints = static_cast<int>(covariates.size()/Nobservables);
}










// MultiLinearModel::MultiLinearModel()
//
// PURPOSE: 
//      Destructor.
//

MultiLinearModel::~MultiLinearModel()
{

}










// MultiLinearModel::getNobservables()
//
// PURPOSE: 
//      Get the protected data member Nobservables.
//
// OUTPUT:
//      An integer containing the number of independent
//      covariates used in the model.
//

int MultiLinearModel::getNobservables()
{
    return Nobservables;
}










// MultiLinearModel::getNpoints()
//
// PURPOSE: 
//      Get the protected data member Npoints.
//
// OUTPUT:
//      An integer containing the number of data points.
//

int MultiLinearModel::getNpoints()
{
    return Npoints;
}












// MultiLinearModel::predict()
//
// PURPOSE:
//      Builds the predictions from a multilinear model of the type y = offset + a*x_1 + b*x_2 + c*x_3 ...
//      where offset, a, b, c, ... are the free parameters and x_i the independent covariates.
//      It is implicit that the covariates are expressed as natural logarithms, otherwise the
//      multilinear model is not valid. The order of the free parameters has to respect the same
//      column order as the input dataset. The offset parameter is always the last one.
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

void MultiLinearModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    // Compute predictions using a loop over the different covariates

    for (int observable = 0; observable < Nobservables; ++observable)
    {
        predictions += covariates.segment(observable*Npoints, Npoints)*modelParameters(observable);
    }


    // Add natural logarithm of proportionality term (offset)
    
    predictions += modelParameters(Nobservables);
}

    






