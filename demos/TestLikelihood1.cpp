#include "TestLikelihood1.h"


// TestLikelihood1::TestLikelihood1()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
// 

TestLikelihood1::TestLikelihood1(const RefArrayXd observations, Model &model)
: Likelihood(observations, model)
{
}









// TestLikelihood1::~TestLikelihood1()
//
// PURPOSE: 
//      Derived class destructor.
//

TestLikelihood1::~TestLikelihood1()
{
}








// TestLikelihood1::logValue()
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

double TestLikelihood1::logValue(RefArrayXd nestedSampleOfParameters)
{
    ArrayXd xCentroid(5);
    xCentroid << -0.400, -0.350, -0.200, 0.100, 0.450;
    ArrayXd yCentroid(5);
    yCentroid << -0.400, 0.200, 0.150, -0.150, 0.100;
    ArrayXd sigma(5);
    sigma << 0.010, 0.010, 0.030, 0.020, 0.050;
    ArrayXd amplitude(5);
    amplitude << 0.500, 1.000, 0.800, 0.500, 0.600;

    ArrayXd logGaussian(5);
    double component1;
    double component2;

    component1 = Functions::logGaussProfile(nestedSampleOfParameters(0), xCentroid(0), sigma(0), amplitude(0));
    component2 = Functions::logGaussProfile(nestedSampleOfParameters(1), yCentroid(0), sigma(0), amplitude(0));
    logGaussian(0) = component1 + component2;
    double logTotalLikelihood = logGaussian(0);
    
    for (int i=1; i < 5; i++)
    {
        component1 = Functions::logGaussProfile(nestedSampleOfParameters(0), xCentroid(i), sigma(i), amplitude(i));
        component2 = Functions::logGaussProfile(nestedSampleOfParameters(1), yCentroid(i), sigma(i), amplitude(i));
        logGaussian(i) = component1 + component2;
        logTotalLikelihood = Functions::logExpSum(logTotalLikelihood,logGaussian(i));
    }
    

    return logTotalLikelihood;
} // END TestLikelihood1::logValue()



