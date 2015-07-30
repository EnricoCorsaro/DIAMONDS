// Derived class for mean likelihood computations
// Created by Enrico Corsaro @ IvS - 19 February 2013
// e-mail: emncorsaro@gmail.com
// Header file "MeanNormalLikelihood.h"
// Implementations contained in "MeanNormalLikelihood.cpp"


#ifndef MEANNORMALLIKELIHOOD_H
#define MEANNORMALLIKELIHOOD_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include "Likelihood.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class MeanNormalLikelihood : public Likelihood
{

    public:

        MeanNormalLikelihood(const RefArrayXd observations, const RefArrayXd uncertainties, Model &model);
        ~MeanNormalLikelihood();
        ArrayXd getUncertainties();
        ArrayXd getNormalizedUncertainties();
        ArrayXd getWeights();

        virtual double logValue(RefArrayXd const modelParameters);

    private:
        
        ArrayXd uncertainties;
        ArrayXd normalizedUncertainties;
        ArrayXd weights;

}; // END class MeanNormalLikelihood

#endif
