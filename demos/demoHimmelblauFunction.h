
#ifndef HIMMELBLAULIKELIHOOD_H
#define HIMMELBLAULIKELIHOOD_H

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

class HimmelblauLikelihood : public Likelihood
{

    public:

        HimmelblauLikelihood(const RefArrayXd observations, Model &model);
        ~HimmelblauLikelihood();

        virtual double logValue(RefArrayXd nestedSampleOfParameters);


    private:

}; 




// HimmelblauLikelihood::HimmelblauLikelihood()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
// 

HimmelblauLikelihood::HimmelblauLikelihood(const RefArrayXd observations, Model &model)
: Likelihood(observations, model)
{
}









// HimmelblauLikelihood::~HimmelblauLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

HimmelblauLikelihood::~HimmelblauLikelihood()
{
}








// HimmelblauLikelihood::logValue()
//
// PURPOSE:
//      Compute a Toy Likelihood example: Himmelblau's function
//
// INPUT:
//      nestedSampleOfParameters: a one-dimensional array containing the actual
//                                values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the log-likelihood value.
//

double HimmelblauLikelihood::logValue(RefArrayXd nestedSampleOfParameters)
{
    double x = nestedSampleOfParameters(0);
    double y = nestedSampleOfParameters(1);
    double logTotalLikelihood = -1.0*pow((x*x + y - 11.0),2) -1.0*pow((x + y*y - 7.0),2);
    
    return logTotalLikelihood;
}


#endif
