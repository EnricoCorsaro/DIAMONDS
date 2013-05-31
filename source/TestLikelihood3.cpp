#include "TestLikelihood3.h"


// TestLikelihood3::TestLikelihood3()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
// 

TestLikelihood3::TestLikelihood3(const RefArrayXd observations, const RefArrayXd uncertainties, Model &model)
: Likelihood(observations, uncertainties, model)
{
} // END TestLikelihood3::TestLikelihood3()









// TestLikelihood3::~TestLikelihood3()
//
// PURPOSE: 
//      Derived class destructor.
//

TestLikelihood3::~TestLikelihood3()
{

} // END TestLikelihood3::~TestLikelihood3()








// TestLikelihood3::logValue()
//
// PURPOSE:
//      Compute the Toy Likelihood example #2 by Feroz et al. 2008, MNRAS, 384, 449.
//
// INPUT:
//      nestedSampleOfParameters: a one-dimensional array containing the actual
//      values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the log-likelihood value.
//

double TestLikelihood3::logValue(RefArrayXd nestedSampleOfParameters)
{
    double logLikelihood;
    double xCenter1 = -3.5;
    double yCenter1 = 0.0;
    double xCenter2 = 3.5;
    double yCenter2 = 0.0;
    double sigma1 = 0.1;
    double sigma2 = 0.1;
    double radius1 = 2.0;
    double radius2 = 2.0;
    double component1;
    double component2;

    component1 = -0.5 * log(2*Functions::PI*sigma1*sigma1) -0.5* pow(sqrt( pow(nestedSampleOfParameters(0) - xCenter1,2) + 
    pow(nestedSampleOfParameters(1) - yCenter1,2) ) - radius1,2)/(sigma1*sigma1);
    component2 = -0.5 * log(2*Functions::PI*sigma2*sigma2) -0.5* pow(sqrt( pow(nestedSampleOfParameters(0) - xCenter2,2) + 
    pow(nestedSampleOfParameters(1) - yCenter2,2) ) - radius2,2)/(sigma2*sigma2);

    logLikelihood = Functions::logExpSum(component1,component2); 

    return logLikelihood;
}



