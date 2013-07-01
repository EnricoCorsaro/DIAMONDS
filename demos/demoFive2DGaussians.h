
#ifndef FIVE2DGAUSSIANS_H
#define FIVE2DGAUSSIANS_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "Functions.h"
#include <Eigen/Core>
#include "Likelihood.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


// Definition of the class

class Multiple2DGaussiansLikelihood : public Likelihood
{
    public:

        Multiple2DGaussiansLikelihood(const RefArrayXd observations, Model &model);
        ~Multiple2DGaussiansLikelihood();

        virtual double logValue(RefArrayXd modelParameters);


    private:
}; 







// Multiple2DGaussiansLikelihood::Multiple2DGaussiansLikelihood()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
// 

Multiple2DGaussiansLikelihood::Multiple2DGaussiansLikelihood(const RefArrayXd observations, Model &model)
: Likelihood(observations, model)
{
}









// Multiple2DGaussiansLikelihood::~Multiple2DGaussiansLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

Multiple2DGaussiansLikelihood::~Multiple2DGaussiansLikelihood()
{
}








// Multiple2DGaussiansLikelihood::logValue()
//
// PURPOSE:
//      Compute the Toy Likelihood example #1 by Feroz et al. 2008, MNRAS, 384, 449.
//
// INPUT:
//      nestedSampleOfParameters: a one-dimensional array containing the actual
//                                values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the log-likelihood value.
//

double Multiple2DGaussiansLikelihood::logValue(RefArrayXd nestedSampleOfParameters)
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
} 


#endif
