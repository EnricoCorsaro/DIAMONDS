
#ifndef ROSENBROCKLIKELIHOOD_H
#define ROSENBROCKLIKELIHOOD_H

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

class RosenbrockLikelihood : public Likelihood
{

    public:

        RosenbrockLikelihood(const RefArrayXd observations, Model &model);
        ~RosenbrockLikelihood();

        virtual double logValue(RefArrayXd nestedSampleOfParameters);


    private:

}; 




// RosenbrockLikelihood::RosenbrockLikelihood()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
// 

RosenbrockLikelihood::RosenbrockLikelihood(const RefArrayXd observations, Model &model)
: Likelihood(observations, model)
{
}









// RosenbrockLikelihood::~RosenbrockLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

RosenbrockLikelihood::~RosenbrockLikelihood()
{
}








// RosenbrockLikelihood::logValue()
//
// PURPOSE:
//      Compute a Toy Likelihood example: Rosenbrock's function
//
// INPUT:
//      nestedSampleOfParameters: a one-dimensional array containing the actual
//                                values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the log-likelihood value.
//

double RosenbrockLikelihood::logValue(RefArrayXd nestedSampleOfParameters)
{
    double theta1 = nestedSampleOfParameters(0);
    double theta2 = nestedSampleOfParameters(1);
    double logTotalLikelihood = -1.0*((1.0-theta1)*(1.0-theta1) + 100.0*(theta2 - theta1*theta1)*(theta2 - theta1*theta1));
    
    return logTotalLikelihood;
}


#endif
