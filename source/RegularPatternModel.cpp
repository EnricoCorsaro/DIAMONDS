#include "RegularPatternModel.h"


// RegularPatternModel::RegularPatternModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates: one-dimensional array containing the values
//      of the independent variable.
//      Norders: the number of expected radial orders to be fitted
//

RegularPatternModel::RegularPatternModel(const RefArrayXd covariates, const int Norders)
: Model(covariates),
  Norders(Norders)
{

}







// RegularPatternModel::RegularPatternModel()
//
// PURPOSE: 
//      Destructor.
//

RegularPatternModel::~RegularPatternModel()
{

}









// RegularPatternModel::buildRegularPatternModel1()
//
// PURPOSE:
//      Builds the predictions from a RegularPattern model.
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
//      (1) Frequency of maximum power excess (nuMax)
//      (2) First large separation (DeltaNu)
//      (3) First small frequency separation 02 (deltaNu02)
//      (4) First small frequency separation 01 (deltaNu01)
//      (5) ...

void RegularPatternModel::predict(RefArrayXd predictions, const RefArrayXd nestedSampleOfParameters)
{
    Nparameters = nestedSampleOfParameters.size();

    double nuMax = nestedSampleOfParameters(0);
    double DeltaNu = nestedSampleOfParameters(1);
    double deltaNu02 = nestedSampleOfParameters(2);
    double deltaNu01 = nestedSampleOfParameters(3);

    /*
    ArrayXd DeltaNu(Norders);
    ArrayXd deltaNu02(Norders);
    ArrayXd deltaNu01(Norders);
    
    for (int i = 0; i < Norders; i++)
    {
        DeltaNu = nestedSampleOfParameters(i);
        deltaNu02 = nestedSampleOfParameters(Norders + i);
        deltaNu01 = nestedSampleOfParameters(2*Norders + i);
    }
    */

    ArrayXd modes0 = ArrayXd::Zero(predictions.size());
    ArrayXd modes0copy = ArrayXd::Zero(predictions.size());
    ArrayXd modes1 = ArrayXd::Zero(predictions.size());
    ArrayXd modes1copy = ArrayXd::Zero(predictions.size());
    ArrayXd modes2 = ArrayXd::Zero(predictions.size());
    ArrayXd modes2copy = ArrayXd::Zero(predictions.size());
    ArrayXd centroid0(Norders);
    ArrayXd centroid1(Norders);
    ArrayXd centroid2(Norders);
    double height = 16.0;
    double linewidth = 0.4;

    for (int i = 0; i < Norders; i++)
    {
        centroid0[i] = nuMax -(Norders/2.)*DeltaNu + i*DeltaNu;
        centroid1[i] = centroid0[i] + i*DeltaNu/2. - deltaNu01;
        centroid2[i] = centroid0[i] + i*DeltaNu - deltaNu02;
        Functions::modeProfile(modes0copy, covariates, centroid0[i], height, linewidth);
        Functions::modeProfile(modes1copy, covariates, centroid1[i], height, linewidth);
        Functions::modeProfile(modes2copy, covariates, centroid2[i], height, linewidth);
        modes0 += modes0copy;
        modes1 += modes1copy;
        modes2 += modes2copy;
    }

    predictions = modes0 + modes1 + modes2;
}







// RegularPatternModel::getNparameters()
//
// PURPOSE: 
//      Get the private data member Nparameters;
//
// OUTPUT:
//      Returns an integer containing the total number of 
//      free parameters used in the model.
//

int RegularPatternModel::getNparameters()
{
    return Nparameters;
}
