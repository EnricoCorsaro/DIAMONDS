#include "SuperGaussianModel.h"


// SuperGaussianModel::SuperGaussianModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:        one-dimensional array containing the values
//                         of the independent variable of the input data file.
//

SuperGaussianModel::SuperGaussianModel(const RefArrayXd covariates)
: Model(covariates)
{
}










// SuperGaussianModel::SuperGaussianModel()
//
// PURPOSE: 
//      Destructor.
//

SuperGaussianModel::~SuperGaussianModel()
{

}










// SuperGaussianModel::predict()
//
// PURPOSE:
//      Builds the predictions from a model incorporating a super Gaussian function with free parameters
//      mean, half width of the plataeu region, standard deviation of the tails, and amplitude.
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

void SuperGaussianModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    // Define useful variables;

    ArrayXd position;
    position.resize(predictions.size());
    position.setZero();
    double argument; 

    
    // The free parameters are sorted below

    double meanValue = modelParameters(0);
    double halfWidthOfPlateau = modelParameters(1);
    double standardDeviation = modelParameters(2);
    double amplitude = modelParameters(3);
    
    position = covariates - meanValue;

    // Perform a loop over all the covariates to compute proper value of the super Gaussian function

    for (int i=0; i < covariates.size(); ++i)
    {
        if (fabs(position(i)) < halfWidthOfPlateau)
        {   
            // If the point is inside the region of the plateau, set its value to the exact amplitude of the super Gaussian
            
            predictions(i) = amplitude;      
        }
        else
        {
            // If the point is outside the region of the plateau, set its prediction value to that of a Gaussian function 
            
            argument = -0.5 * pow((fabs(position(i)) - halfWidthOfPlateau)/standardDeviation,2);
            predictions(i) = exp(argument)*amplitude;
        }
    }
}
