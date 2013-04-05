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
#include "Functions.h"
#include "Prior.h"
#include "Likelihood.h"


using Eigen::ArrayXd;
using Eigen::ArrayXXd;


class NestedSampler
{

    public:
        
        ArrayXXd posteriorSample;                // parameter values sampling the posterior
        ArrayXd logLikelihoodOfPosteriorSample;  // logLikelihood values corresponding to the posterior sample 
        ArrayXd logWeightOfPosteriorSample;      // logWeights corresponding to the posterior sample

        NestedSampler(Prior &prior, Likelihood &likelihood); 
        ~NestedSampler();
        double getLogEvidence();
        double getLogEvidenceError();
        double getInformationGain();
        int getNiterations();
        void run(const int Nobjects);


	private:

        mt19937 engine;
        double informationGain;
        double logEvidence;
        double logEvidenceError;
        static unsigned nestedCounter;           // Static counter containing the number of Nested processes running
        int Niterations;                         // Counter saving the number of nested loops used
        ArrayXXd nestedSampleOfParameters;       // parameters values (the free parameters of the problem)
        ArrayXd logLikelihood;                   // log-likelihood values corresponding to parameter values
        Prior &prior;
        Likelihood &likelihood;


}; // END class NestedSampler

#endif
