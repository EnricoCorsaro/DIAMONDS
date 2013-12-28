
#ifndef RASTRIGINLIKELIHOOD_H
#define RASTRIGINLIKELIHOOD_H

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

class RastriginLikelihood : public Likelihood
{

    public:

        RastriginLikelihood(const RefArrayXd observations, Model &model);
        ~RastriginLikelihood();

        virtual double logValue(RefArrayXd nestedSampleOfParameters);


    private:

}; 




// RastriginLikelihood::RastriginLikelihood()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
// 

RastriginLikelihood::RastriginLikelihood(const RefArrayXd observations, Model &model)
: Likelihood(observations, model)
{
}









// RastriginLikelihood::~RastriginLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

RastriginLikelihood::~RastriginLikelihood()
{
}








// RastriginLikelihood::logValue()
//
// PURPOSE:
//      Compute a Toy Likelihood example: Rastrigin's function
//
// INPUT:
//      nestedSampleOfParameters: a one-dimensional array containing the actual
//                                values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the log-likelihood value.
//

double RastriginLikelihood::logValue(RefArrayXd nestedSampleOfParameters)
{
    double theta1 = nestedSampleOfParameters(0);
    double theta2 = nestedSampleOfParameters(1);
    double term1 = theta1*theta1 -10.0*cos(2.0*Functions::PI*theta1);
    double term2 = theta2*theta2 -10.0*cos(2.0*Functions::PI*theta2);
    double logTotalLikelihood = -1.0*(10*nestedSampleOfParameters.size() + term1 + term2);
    
    return logTotalLikelihood;
}


#endif
