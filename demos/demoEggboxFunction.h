
#ifndef EGGBOXLIKELIHOOD_H
#define EGGBOXLIKELIHOOD_H

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

class EggboxLikelihood : public Likelihood
{

    public:

        EggboxLikelihood(const RefArrayXd observations, Model &model);
        ~EggboxLikelihood();

        virtual double logValue(RefArrayXd modelParameters);


    private:

};





// EggboxLikelihood::EggboxLikelihood()
//
// PURPOSE: 
//      Derived class onstructor.
//
// INPUT:
// 

EggboxLikelihood::EggboxLikelihood(const RefArrayXd observations, Model &model)
: Likelihood(observations, model)
{
}









// EggboxLikelihood::~EggboxLikelihood()
//
// PURPOSE: 
//      Derived class destructor.
//

EggboxLikelihood::~EggboxLikelihood()
{
}








// EggboxLikelihood::logValue()
//
// PURPOSE:
//      Compute the Toy Likelihood example #2 by Feroz et al. 2008, MNRAS, 384, 449.
//
// INPUT:
//      nestedSampleOfParameters: a one-dimensional array containing the actual
//                                values of the free parameters that describe the model.
//
// OUTPUT:
//      a double number containing the log-likelihood value.
//

double EggboxLikelihood::logValue(RefArrayXd nestedSampleOfParameters)
{
    double logLikelihood;
    logLikelihood = pow(2.0 + cos(nestedSampleOfParameters(0)/2.0) * cos(nestedSampleOfParameters(1)/2.0), 5); 

    return logLikelihood;
}


#endif
