// Class for nested sampling inference
// Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "NestedSampler.h"
// Implementation contained in "NestedSampler.cpp"

#ifndef NESTEDSAMPLER_H
#define NESTEDSAMPLER_H

#include <cfloat>
#include <random>
#include <Eigen/Dense>
#include "MathExtra.h"
#include "RandomVariate.h"

using Eigen::ArrayXd;
using Eigen::ArrayXXd;


class NestedSampler
{
    public:
        
        NestedSampler(RandomVariate &variate); 
        ~NestedSampler();
        void run(int Nobjects, int Niter);
        double getLogEvidence();
        double getLogEvidenceError();
        double getInformationH();
        ArrayXd parameterObjects;               // parameter values (the free parameters of the problem)
        ArrayXd posteriorSample;                // parameter values sampled from the posterior
        ArrayXd logLikelihoodOfPosteriorSample; // logLikelihood values corresponding to the posterior sample 

	private:
        
        double informationH;
        double logEvidence;
        double logEvidenceError;
        ArrayXd logLikelihood;                  // log-likelihood values corresponding to parameter values
        ArrayXd logWeight;                      // sum(weight) = Evidence Z
        static long nestedCounter;              // Static index containing the number of Nested processes running
}; // END class NestedSampler

#endif
