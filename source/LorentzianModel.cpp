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
//      nestedSampleOfParameters: one-dimensional array where each element
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

void LorentzianModel::predict(RefArrayXd predictions, const RefArrayXd nestedSampleOfParameters)
{
    Nparameters = nestedSampleOfParameters.size();

    if (Nparameters == 1)
    {
        Functions::lorentzProfile(predictions, covariates, nestedSampleOfParameters(0));
    }
    else 
        if (Nparameters == 2)
        {
            Functions::lorentzProfile(predictions, covariates, nestedSampleOfParameters(0), nestedSampleOfParameters(1));
        }
    else 
        if (Nparameters == 3)
        {
            Functions::lorentzProfile(predictions, covariates, nestedSampleOfParameters(0), nestedSampleOfParameters(1), nestedSampleOfParameters(2));
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
