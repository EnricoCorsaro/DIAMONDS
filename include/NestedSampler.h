// Class for nested sampling inference
// Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "NestedSampler.h"
// Implementation contained in "NestedSampler.cpp"

#ifndef NESTEDSAMPLER_H
#define NESTEDSAMPLER_H

#include <iostream>
#include <iomanip>
#include <cfloat>
#include <random>
#include <ctime>
#include <cmath>
#include <vector>
#include <cassert>
#include <Eigen/Dense>
#include "Functions.h"
#include "Prior.h"
#include "Likelihood.h"
#include "Metric.h"
#include "Clusterer.h"


using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXi> RefArrayXi;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;

class NestedSampler
{

    public:
        
        ArrayXXd posteriorSample;                // parameter values sampling the posterior
        ArrayXd logLikelihoodOfPosteriorSample;  // logLikelihood values corresponding to the posterior sample 
        ArrayXd logWeightOfPosteriorSample;      // logWeights corresponding to the posterior sample

        NestedSampler(Prior &prior, Likelihood &likelihood, Metric &metric, Clusterer &clusterer); 
        ~NestedSampler();
        
        double getLogEvidence();
        double getLogEvidenceError();
        double getInformationGain();
        int getNiterations();
        void run(const int Nobjects, const int NiterationsBeforeClustering, const int Ndraws = 1);
        virtual void drawWithConstraint(const RefArrayXXd totalSampleOfParameters, const int Nclusters, const RefArrayXi clusterIndices,
                                        const double logWidthInPriorMass, RefArrayXXd drawnSampleOfParameters) = 0;


    protected:

        Prior &prior;
        Likelihood &likelihood;
        Metric &metric;
        Clusterer &clusterer;
        mt19937 engine;
        int Ndimensions;


	private:

        double informationGain;
        double logEvidence;
        double logEvidenceError;
        static unsigned nestedCounter;           // Static counter containing the number of Nested processes running
        int Niterations;                         // Counter saving the number of nested loops used
        ArrayXXd nestedSampleOfParameters;       // parameters values (the free parameters of the problem)
        ArrayXd logLikelihood;                   // log-likelihood values corresponding to parameter values


}; // END class NestedSampler

#endif
