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
#include "Prior.h"
#include "Likelihood.h"


using Eigen::ArrayXd;
using Eigen::ArrayXXd;


class NestedSampler
{

    public:
        
        ArrayXXd posteriorSample;                // parameter values sampled from the posterior
        ArrayXd logLikelihoodOfPosteriorSample;  // logLikelihood values corresponding to the posterior sample 
        
        NestedSampler(Prior &prior, Likelihood &likelihood); 
        ~NestedSampler();
        double getLogEvidence();
        double getLogEvidenceError();
        double getInformationGain();
        int getNestIteration();
        void run();


	private:
        
        double informationGain;
        double logEvidence;
        double logEvidenceError;
        static unsigned nestedCounter;           // Static counter containing the number of Nested processes running
        int nestIteration;                       // Counter saving the number of nested loops used
        ArrayXXd nestedParameters;               // parameters values (the free parameters of the problem)
        ArrayXd logLikelihood;                   // log-likelihood values corresponding to parameter values
        ArrayXd logWeight;                       // log(prior mass * Likelihood), accumulating evidence
        Prior &prior;
        Likelihood &likelihood;


}; // END class NestedSampler

#endif
