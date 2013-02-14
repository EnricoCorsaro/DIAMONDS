// Class for computing likelihoods given
// the model and the data.
// Enrico Corsaro @ IvS - 14 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "Likelihood.h"
// Implementation contained in "Likelihood.cpp"

#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <Eigen/Core>
#include "MathExtra.h"
#include "Model.h"
#include "File.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;

class Likelihood
{

    public:

        Likelihood(const RefArrayXd data, const int modelID);
        void setLikelihood(const RefArrayXd modelParameters)
        double buildLogGaussianLikelihood();
        double buildMarginalizedLogGaussianLikelihood();
        double buildMedianLikelihood();
        ArrayXd getCovariates();
        ArrayXd getObservations();
        ArrayXd getUncertainties();
        ArrayXd getPredictions();


    private:

        int modelIdentifier;
        ArrayXd covariates;
        ArrayXd observations;
        ArrayXd uncertainties;
        ArrayXd predictions;

}; // END class Likelihood

#endif
