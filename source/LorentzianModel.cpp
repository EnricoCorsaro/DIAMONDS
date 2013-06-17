#include "LorentzianModel.h"


// LorentzianModel::LorentzianModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates: one-dimensional array containing the values
//      of the independent variable.
//

LorentzianModel::LorentzianModel(const RefArrayXd covariates)
: Model(covariates)
{

}







// LorentzianModel::LorentzianModel()
//
// PURPOSE: 
//      Destructor.
//

LorentzianModel::~LorentzianModel()
{

}









// LorentzianModel::buildLorentzianModel1()
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
//      (1) centroid of the Lorentz profile
//      (2) amplitude of the Lorentz profile (defaulted to 1.0 if not given)
//      (3) width of the Lorentz profile (defaulted to 1.0 if not given)
//

void LorentzianModel::predict(RefArrayXd predictions, const RefArrayXd modelParameters)
{
    Nparameters = modelParameters.size();

    if (Nparameters == 1)
    {
        Functions::lorentzProfile(predictions, covariates, modelParameters(0));
    }
    else 
        if (Nparameters == 2)
        {
            Functions::lorentzProfile(predictions, covariates, modelParameters(0), modelParameters(1));
        }
    else 
        if (Nparameters == 3)
        {
            Functions::lorentzProfile(predictions, covariates, modelParameters(0), modelParameters(1), modelParameters(2));
        }
    else
    {
        cerr << "Number of free parameters does not match model Lorentzian. Quitting program." << endl;
        exit(EXIT_FAILURE);
    }
}







// LorentzianModel::getNparameters()
//
// PURPOSE: 
//      Get the private data member Nparameters;
//
// OUTPUT:
//      Returns an integer containing the total number of 
//      free parameters used in the model.
//

int LorentzianModel::getNparameters()
{
    return Nparameters;
}
