#include "TestLikelihood2.h"


// TestLikelihood2::TestLikelihood2()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
// 

TestLikelihood2::TestLikelihood2(const RefArrayXd observations, const RefArrayXd uncertainties, Model &model)
: Likelihood(observations, uncertainties, model)
{
} // END TestLikelihood2::TestLikelihood2()









// TestLikelihood2::~TestLikelihood2()
//
// PURPOSE: 
//      Derived class destructor.
//

TestLikelihood2::~TestLikelihood2()
{

} // END TestLikelihood2::~TestLikelihood2()








// TestLikelihood2::logValue()
//
// PURPOSE:
//      Compute the Toy Likelihood example #1 by Feroz et al. 2008, MNRAS, 384, 449.
//
// INPUT:
//      nestedSampleOfParameters: a one-dimensional array containing the actual
//      values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the log-likelihood value.
//

double TestLikelihood2::logValue(RefArrayXd nestedSampleOfParameters)
{
    double logLikelihood;
    logLikelihood = pow(2.0 + cos(nestedSampleOfParameters(0)/2.0) * cos(nestedSampleOfParameters(1)/2.0), 5); 

    return logLikelihood;
}



