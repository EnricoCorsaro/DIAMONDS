
#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <Eigen/Core>
#include "MathExtra.h"
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;



class Likelihood
{

    public:

        Likelihood(RefArrayXd covariates, RefArrayXd observations, RefArrayXd uncertainties, Model &model);
        virtual double logDensity(RefArrayXd modelParameters) = 0;

    protected:

    private:

        ArrayXd covariates;
        ArrayXd observations;
        ArrayXd uncertainties;
        Model &model;
        
};

#endif
