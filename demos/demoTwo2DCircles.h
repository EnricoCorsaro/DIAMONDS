
#ifndef DEMOTWOCIRCLES_H
#define DEMOTWOCIRCLES_H

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


// Definition of the Likelihood

class TwoCirclesLikelihood : public Likelihood
{

    public:

        TwoCirclesLikelihood(const RefArrayXd observations, Model &model);
        ~TwoCirclesLikelihood();

        virtual double logValue(RefArrayXd modelParameters);


    private:

}; // END class TwoCirclesLikelihood






// TwoCirclesLikelihood::TwoCirclesLikelihood()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
// 

TwoCirclesLikelihood::TwoCirclesLikelihood(const RefArrayXd observations, Model &model)
: Likelihood(observations, model)
{
}









// TwoCirclesLikelihood::~TwoCirclesLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

TwoCirclesLikelihood::~TwoCirclesLikelihood()
{
}








// TwoCirclesLikelihood::logValue()
//
// PURPOSE:
//      Compute the Toy Likelihood example #3 by Feroz et al. 2008, MNRAS, 384, 449.
//
// INPUT:
//      nestedSampleOfParameters: a one-dimensional array containing the actual
//                                values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the log-likelihood value.
//

double TwoCirclesLikelihood::logValue(RefArrayXd nestedSampleOfParameters)
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




#endif
