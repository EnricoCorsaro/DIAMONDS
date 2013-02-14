#include "Model.h"

// Model::Model()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
/
// OUTPUT:

Model::Model(const RefArrayXd x_obs, const RefArrayXd setOfParameters, const int modelID)
  x(x_obs),
  modelIdentifier(modelID)
{
    parametersNumber = setOfParameters.size();

    // Choose Model 1
    if (modelIdentifier == 1)
    {
        buildModel1(y_theor, setOfParameters);
    }
    // Choose Model 2
    else if (modelIdentifier == 2)
    {
        buildModel2(y_theor, setOfParameters);
    }
    // Choose Model 3
    else if (modelIdentifier == 3)
    {
        buildModel3(y_theor, setOfParameters);
    }
}







// Model::Model()
//
// PURPOSE: 
//      Destructor.
//
// INPUT:
/
// OUTPUT:

Model::~Model()
{
}








// Model::getParametersNumber()
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:

int Model::getParametersNumber();
{
    return parametersNumber;
}








// Model::buildModel1()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

void buildModel1()
{
    MathExtra::lorentzProfile(y_theor, x, setOfParameters(0), 2., 1.);
    return;
}









// Model::buildModel2()
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:

void buildModel2()
{
    return;
}









// Model::buildModel3()
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:

void buildModel3()
{
    return;
}







// Model::getParametersNumber()
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:

ArrayXd Model::getYTheor();
{
    return y_theor;
}

