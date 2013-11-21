
#ifndef SINGLENDGAUSSIANLIKELIHOOD_H
#define SINGLENDGAUSSIANLIKELIHOOD_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <functional>
#include "Functions.h"
#include <Eigen/Core>
#include "Likelihood.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


// Class definition

class SingleNDGaussianLikelihood : public Likelihood
{

    public:

        SingleNDGaussianLikelihood(const RefArrayXd observations, Model &model,int Ndimensions);
        ~SingleNDGaussianLikelihood();

        virtual double logValue(RefArrayXd nestedSampleOfParameters);


    private:
    
        int Ndimensions;
        ArrayXd centroid;
        ArrayXd sigma;

}; 




// SingleNDGaussianLikelihood::SingleNDGaussianLikelihood()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
// 

SingleNDGaussianLikelihood::SingleNDGaussianLikelihood(const RefArrayXd observations, Model &model, int Ndimensions)
: Likelihood(observations, model),
  Ndimensions(Ndimensions)
{
    centroid.resize(Ndimensions);
    sigma.resize(Ndimensions);

    mt19937 engine;
    clock_t clockticks = clock();
    engine.seed(clockticks);

    uniform_real_distribution<double> uniform1(0.1,0.5);
    uniform_real_distribution<double> uniform2(-10.0,10.0);

    cout << "Gaussian Likelihood in " << Ndimensions << " dimensions" << endl;

    for (int i=0; i < Ndimensions; ++i)
    {
        centroid(i) = uniform2(engine);
        sigma(i) = uniform1(engine);
        cout << "Dimension #: " << i << "   " << "Centroid: " << centroid(i) << "   " << "sigma: " << sigma(i) << endl;
    }

}









// SingleNDGaussianLikelihood::~SingleNDGaussianLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

SingleNDGaussianLikelihood::~SingleNDGaussianLikelihood()
{
}








// SingleNDGaussianLikelihood::logValue()
//
// PURPOSE:
//      Compute a Toy Likelihood example: Ndiemsional Gaussian
//
// INPUT:
//      nestedSampleOfParameters: a one-dimensional array containing the actual
//                                values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the log-likelihood value.
//

double SingleNDGaussianLikelihood::logValue(RefArrayXd nestedSampleOfParameters)
{
    assert(nestedSampleOfParameters.size() == Ndimensions);
    ArrayXd exponent(Ndimensions);
   
    exponent = ((nestedSampleOfParameters - centroid)/sigma).square();
    
    double logTotalLikelihood = -log(2*Functions::PI*sigma.prod()) -0.5*exponent.sum();
    return logTotalLikelihood;
}


#endif
