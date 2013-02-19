#include "SingleLorentzianModel.h"


// SingleLorentzianModel::SingleLorentzianModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates: one-dimensional array containing the values
//      of the independent variable.
//

SingleLorentzianModel::SingleLorentzianModel(const RefArrayXd covariates)
: Model(covariates)
{

} // END SingleLorentzianModel::SingleLorentzianModel








// SingleLorentzianModel::SingleLorentzianModel()
//
// PURPOSE: 
//      Destructor.
//

SingleLorentzianModel::~SingleLorentzianModel()
{

} // END SingleLorentzianModel::~SingleLorentzianModel










// SingleLorentzianModel::buildSingleLorentzianModel1()
//
// PURPOSE:
//      Builds the predictions from a single Lorentzian model.
//
// INPUT:
//      predictions: one-dimensional array to contain the predictions
//      from the model
//      modelParameters: one-dimensional array where each element
//      contains the value of one free parameter of the model
//
// OUTPUT:
//      void
//
// NOTE:
//      The free parameters are to be given in the order
//      (1) centroid of the peak (frequency)
//      (2) amplitude of the peak (observed flux)
//      (3) mode linewidth (related to the lifetime of the mode)
//

void SingleLorentzianModel::predict(RefArrayXd predictions, const RefArrayXd modelParameters)
{
    if (modelParameters.size() == 1)
    {
    MathExtra::lorentzProfile(predictions, covariates, modelParameters(0));
    }
    else if (modelParameters.size() == 2)
    {
    MathExtra::lorentzProfile(predictions, covariates, modelParameters(0), modelParameters(1));
    }
    else if (modelParameters.size() == 3)
    {
    MathExtra::lorentzProfile(predictions, covariates, modelParameters(0), modelParameters(1), modelParameters(2));
    }
    else
    {
        cerr "Number of free parameters do not match model Single Lorentzian. Quitting program." << endl;
        exit(1);
    }
} // END SingleLorentzianModel::predict()








// SingleLorentzianModel::getParametersNumber()
//
// PURPOSE: 
//      Get the private data member parametersNumber;
//
// OUTPUT:
//      Returns an integer containing the total number of 
//      free parameters used in the model.
//

int SingleLorentzianModel::getParametersNumber();
{
    return parametersNumber;
} // END SingleLorentzianModel::getParametersNumber()
