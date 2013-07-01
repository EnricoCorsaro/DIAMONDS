
#ifndef DEMOTWO2DGAUSSIANS_H
#define DEMOTWO2DGAUSSIANS_H

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
//      Compute a toy model of two well separated Gaussians. 
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
    ArrayXd xCentroid(2);
    xCentroid << -3.0, +3.0;
    ArrayXd yCentroid(2);
    yCentroid << 0.0, 0.0;
    ArrayXd sigma(2);
    sigma << 1.0, 1.0; 
    ArrayXd amplitude(2);
    amplitude << 0.5, 0.5; 

    ArrayXd logGaussian(2);
 
    for (int i=0; i < 2; i++)
    {
        double component1 = Functions::logGaussProfile(nestedSampleOfParameters(0), xCentroid(i), sigma(i), amplitude(i));
        double component2 = Functions::logGaussProfile(nestedSampleOfParameters(1), yCentroid(i), sigma(i), amplitude(i));
        logGaussian(i) = component1 + component2;
    }
    
    double logTotalLikelihood = Functions::logExpSum(logGaussian(0), logGaussian(1));
    return logTotalLikelihood;
} 


#endif
