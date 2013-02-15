#ifndef NORMALLIKELIHOOD_H
#define NORMALLIKELIHOOD_H

#include <cmath>
#include <Eigen/Core>
#include "MathExtra.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class NormalLikelihood : public Likelihood
{

    public:

        NormalLikelihood(RefArrayXd covariates, RefArrayXd observations, RefArrayXd uncertainties, Model &model);
        ~NormalLikelihood();

        double logDensity(RefArrayXd modelParameters);

    protected:
    
        

    private:

};

#endif
