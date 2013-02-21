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

} // END LorentzianModel::LorentzianModel








// LorentzianModel::LorentzianModel()
//
// PURPOSE: 
//      Destructor.
//

LorentzianModel::~LorentzianModel()
{

} // END LorentzianModel::~LorentzianModel










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
//      (1) centroid of the peak (frequency)
//      (2) amplitude of the peak (observed flux)
//      (3) mode linewidth (related to the lifetime of the mode)
//

void LorentzianModel::predict(RefArrayXd predictions, const RefArrayXd modelParameters)
{
    parametersNumber = modelParameters.size();

    if (parametersNumber == 1)
    {
        MathExtra::lorentzProfile(predictions, covariates, modelParameters(0));
    }
    else 
        if (parametersNumber == 2)
        {
            MathExtra::lorentzProfile(predictions, covariates, modelParameters(0), modelParameters(1));
        }
    else 
        if (parametersNumber == 3)
        {
            MathExtra::lorentzProfile(predictions, covariates, modelParameters(0), modelParameters(1), modelParameters(2));
        }
    else
    {
        cerr << "Number of free parameters does not match model Lorentzian. Quitting program." << endl;
        exit(EXIT_FAILURE);
    }
} // END LorentzianModel::predict()








// LorentzianModel::getParametersNumber()
//
// PURPOSE: 
//      Get the private data member parametersNumber;
//
// OUTPUT:
//      Returns an integer containing the total number of 
//      free parameters used in the model.
//

int LorentzianModel::getParametersNumber()
{
    return parametersNumber;
} // END LorentzianModel::getParametersNumber()
