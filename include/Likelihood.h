// Abstract base class for likelihood computations.
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 15 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "Likelihood.h"
// Implementations contained in "Likelihood.cpp"

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

        Likelihood(const RefArrayXd covariates, const RefArrayXd observations, const RefArrayXd uncertainties, Model &model);
        ~Likelihood();
        ArrayXd getCovariates();
        ArrayXd getObservations();

        virtual double logValue(RefArrayXd nestedSampleOfParameters) = 0;

    protected:
        
        ArrayXd covariates;
        ArrayXd observations;
        ArrayXd uncertainties;
        Model &model;

    private:

        
}; // END class Likelihood

#endif
