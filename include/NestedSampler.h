// Class for nested sampling inference
// Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster,kuleuven.be
// Header file "NestedSampler.h"

#ifndef NESTEDSAMPLER_H
#define NESTEDSAMPLER_H

#include <cfloat>
#include <random>
#include <Eigen/Dense>
#include "MathExtra.h"
#include "FileProcess.h"
#include "RandomVariate.h"

using Eigen::ArrayXd;
using Eigen::ArrayXXd;



class NestedSampler
{
    public:
        
        NestedSampler(RandomVariate &variate); 
        void run(int Nobjects, int Niter);
        double getLogEvidence();
        double getLogEvidenceError();
        double getInformationH();
        ArrayXXd parameter;                     // parameter values (the free parameters of the problem)
        ArrayXXd posteriorSample;               // parameter values sampled from the posterior
        ArrayXd logLikelihoodOfPosteriorSample; // logLikelihood values corresponding to the posterior sample 

	private:
        
        RandomVariate &randomVariate;           // random variate to draw (constrained) prior from
        double informationH;
        double logEvidence;
        double logEvidenceError;
        ArrayXd logLikelihood;                  // log-likelihood values
        ArrayXd logWeight;                      // sum(weight) = Evidence Z
};

#endif
