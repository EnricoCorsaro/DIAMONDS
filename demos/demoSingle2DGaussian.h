
#ifndef SINGLE2DGAUSSIANLIKELIHOOD_H
#define SINGLE2DGAUSSIANLIKELIHOOD_H

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


// Class definition

class Single2DGaussianLikelihood : public Likelihood
{

    public:

        Single2DGaussianLikelihood(const RefArrayXd observations, Model &model);
        ~Single2DGaussianLikelihood();

        virtual double logValue(RefArrayXd nestedSampleOfParameters);


    private:

}; 




// Single2DGaussianLikelihood::Single2DGaussianLikelihood()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
// 

Single2DGaussianLikelihood::Single2DGaussianLikelihood(const RefArrayXd observations, Model &model)
: Likelihood(observations, model)
{
}









// Single2DGaussianLikelihood::~Single2DGaussianLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

Single2DGaussianLikelihood::~Single2DGaussianLikelihood()
{
}








// Single2DGaussianLikelihood::logValue()
//
// PURPOSE:
//      Compute a Toy Likelihood example: 2D Gaussian
//
// INPUT:
//      nestedSampleOfParameters: a one-dimensional array containing the actual
//                                values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the log-likelihood value.
//

double Single2DGaussianLikelihood::logValue(RefArrayXd nestedSampleOfParameters)
{
    double xCentroid = 9.5;
    double yCentroid = 20.0;;
    double xSigma = 1.5;
    double ySigma = 1.5;

    double exponent1 = pow((nestedSampleOfParameters(0) - xCentroid)/xSigma, 2);
    double exponent2 = pow((nestedSampleOfParameters(1) - yCentroid)/ySigma, 2);
    double logTotalLikelihood = -log(2*Functions::PI*xSigma*ySigma) -0.5*(exponent1 + exponent2);
    
    return logTotalLikelihood;
}


#endif