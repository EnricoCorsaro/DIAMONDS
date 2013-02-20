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

} // END Model::Model()








// Model::Model()
//
// PURPOSE: 
//      Destructor.
//

Model::~Model()
{

} // END Model::~Model()








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
} // END Model::getCovariates()


