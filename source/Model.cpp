#include "Model.h"


// Model::Model()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates: one-dimensional array containing the values
//      of the independent variable.
//

Model::Model(const RefArrayXd covariates)
: covariates(covariates)
{

}








// Model::Model()
//
// PURPOSE: 
//      Destructor.
//

Model::~Model()
{

}








// Model::getCovariates()
//
// PURPOSE: 
//      Get the protected data member covariates.
//
// OUTPUT:
//      A one-dimensional array containing the values
//      of the independent variable.
//

ArrayXd Model::getCovariates()
{
    return covariates;
}











// Model::getNparameters()
//
// PURPOSE: 
//      Get the protected data member Nparameters;
//
// OUTPUT:
//      Returns an integer containing the total number of 
//      free parameters used in the model.
//

int Model::getNparameters()
{
    return Nparameters;
}
