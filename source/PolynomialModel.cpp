#include "PolynomialModel.h"


// PolynomialModel::PolynomialModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:                 one-dimensional array containing the values
//                                  of the independent variable.
//      covariatesUncertainties:    ome-dimensional array containing the values of the uncertainties on the independent variable
//      Ndegrees:                   an integer containing the degree of the polynomial to fit.
//                                  A polynomial of the type y(x) = a*x + b has degree = 1.
//      covariatesOffset:           a double specifying the starting point of the intercept in the covariates
//

PolynomialModel::PolynomialModel(const RefArrayXd covariates, const RefArrayXd covariatesUncertainties, const int Ndegrees, const double covariatesOffset)
: Model(covariates),
  covariatesUncertainties(covariatesUncertainties),
  Ndegrees(Ndegrees),
  covariatesOffset(covariatesOffset)
{
    Npoints = covariates.size();
}










// PolynomialModel::PolynomialModel()
//
// PURPOSE: 
//      Destructor.
//

PolynomialModel::~PolynomialModel()
{

}










// PolynomialModel::getNdegrees()
//
// PURPOSE: 
//      Get the protected data member Ndegrees.
//
// OUTPUT:
//      An integer containing the degree of the polynomial.
//

int PolynomialModel::getNdegrees()
{
    return Ndegrees;
}









// PolynomialModel::getCovariatesOffset()
//
// PURPOSE: 
//      Get the protected data member covariatesOffset.
//
// OUTPUT:
//      A double containing the offset used for the covariates.
//

double PolynomialModel::getCovariatesOffset()
{
    return covariatesOffset;
}










// PolynomialModel::predict()
//
// PURPOSE:
//      Builds the predictions from a polynomial model of the type f = offset + a*x + b*x^2 + c*x^3 + ...
//      where offset, a, b, c, ... are the free parameters and x the covariates.
//      The free parameters related to the covariates have to be sorted in increasing degree order.
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

void PolynomialModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    // Compute predictions with a loop over the different degrees of the polynomial

    for (int degree = 0; degree < Ndegrees; ++degree)
    {
        predictions += (covariates - covariatesOffset).pow(degree + 1)*modelParameters(degree);
    }

    // Add offset term to the polynomial
    
    predictions += modelParameters(Ndegrees);
}
















// PolynomialModel::computeVariance()
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

void PolynomialModel::computeVariance(RefArrayXd modelVariance, RefArrayXd const modelParameters)
{
    ArrayXd polynomialSDV(Npoints);
    polynomialSDV.setZero();
    
    // Compute the total uncertainty arising from the polynomial model
    
    for (int degree = 0; degree < Ndegrees; ++degree)
    {
        polynomialSDV += (degree + 1)*modelParameters(degree)*(covariates - covariatesOffset).pow(degree);
    }

    polynomialSDV *= covariatesUncertainties;
    modelVariance += polynomialSDV.square();
} 
    





