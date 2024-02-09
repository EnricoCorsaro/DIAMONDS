#include "MultiLinearModel.h"


// MultiLinearModel::MultiLinearModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:                 one-dimensional array containing the values
//                                  of the independent variables concatenated in the same
//                                  order as the columns from left to right in the input data file.
//      covariatesUncertainties:    ome-dimensional array containing the values of the uncertainties on the independent variables concatenated
//                                  in the same order as the columns from left to right in the input data file
//      Nobservables:               the number of observables in the multilinear model
//

MultiLinearModel::MultiLinearModel(const RefArrayXd covariates, const RefArrayXd covariatesUncertainties, int Nobservables)
: Model(covariates),
  covariatesUncertainties(covariatesUncertainties),
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










// MultiLinearModel::getCovariatesUncertainties();
//
// PURPOSE:
//      Get protected data member covariatesUncertainties.
//
// OUTPUT:
//      covariatesUncertainties: one-dimensional array containing the
//      uncertainties on the independent variables.
//

ArrayXd MultiLinearModel::getCovariatesUncertainties()
{
    return covariatesUncertainties;
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

    







// MultiLinearModel::computeVariance()
//
// PURPOSE:
//      Builds the total covariates uncertainties by applying a standard error propagation from the model equation.
//
// INPUT:
//      modelVariance:        one-dimensional array to contain the total variance on the covariates
//                            as computed by the model
//      modelParameters:      one-dimensional array where each element
//                            contains the value of a free parameter of the model
//
// OUTPUT:
//      void
//

void MultiLinearModel::computeVariance(RefArrayXd modelVariance, RefArrayXd const modelParameters)
{
    // Compute the total uncertainty arising from the multi-linear model 

    for (int observable = 0; observable < Nobservables; ++observable)
    {
        modelVariance += covariatesUncertainties.segment(observable*Npoints, Npoints).square()
                        *modelParameters(observable)*modelParameters(observable);
    }
} 