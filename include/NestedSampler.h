// Class for nested sampling inference
// Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster,kuleuven.be
// Header file "NestedSampler.h"

#ifndef NESTEDSAMPLER_H
#define NESTEDSAMPLER_H

#include <cfloat>
#include "MathExtra.h"
#include "FileProcess.h"
#include "RandomVariate.h"


class NestedSampler
{
    public:
        
        NestedSampler(RandomVariate &variate); 
        void run(int Nobjects, int Niter);
        double getLogEvidence();
        double getLogEvidenceError();
        double getInformationH();
        vector<double> param;                           // parameter values (the free parameters of the problem)
        vector<double> posteriorSample;                 // parameter values sampled from the posterior
        vector<double> logLikelihoodOfPosteriorSample;  // logLikelihood values corresponding to the posterior sample 
        vector<double> results;                         // output logZ, logZ_err, information H

	private:
        
        RandomVariate &randomVariate;           // random variate to draw (constrained) prior from
        double informationH;
        double logEvidence;
        double logEvidenceError;
        vector<double> logLikelihood;           // log-likelihood values
        vector<double> logWeight;               // sum(weight) = Evidence Z
        double updateInformationGain(double H_old, double logZ_old, double logZ_new, int worst);
};

#endif
