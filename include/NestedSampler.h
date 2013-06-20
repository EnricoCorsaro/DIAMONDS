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
        
        ArrayXXd posteriorSample;                // parameter values for sampling the posterior
        ArrayXd logLikelihoodOfPosteriorSample;  // logLikelihood values corresponding to the posterior sample 
        ArrayXd logWeightOfPosteriorSample;      // logWeights corresponding to the posterior sample

        NestedSampler(const bool printOnTheScreen, const int Nobjects, vector<Prior*> ptrPriors, 
                      Likelihood &likelihood, Metric &metric, Clusterer &clusterer); 
        ~NestedSampler();
        
        void run(const double terminationFactor = 0.5, 
                 const int NiterationsBeforeClustering = 10, const int maxNdrawAttempts = 200);
        virtual void drawWithConstraint(const RefArrayXXd totalSampleOfParameters, const int Nclusters, const RefArrayXi clusterIndices,
                                        const double logWidthInPriorMass, RefArrayXXd drawnSampleOfParameters, const int maxNdrawAttempts) = 0;
        int getNiterations();
        double getLogEvidence();
        double getLogEvidenceError();
        double getLogMeanEvidence();
        double getLogMeanEvidenceError();
        double getLogMeanTotalEvidence();
        double getLogMeanTotalEvidenceError();
        double getInformationGain();
        double getComputationalTime();


    protected:

        vector<Prior*> ptrPriors;
        Likelihood &likelihood;
        Metric &metric;
        Clusterer &clusterer;
        bool printOnTheScreen;
        int Ndimensions;
        int Nobjects;                           // Total number of objects
        double actualLogLikelihoodConstraint;   // Likelihood constraining value at each nested iteration
        double logTotalWidthInPriorMass;        // The remaining width in prior mass at a given nested iteration (log X_k)

        mt19937 engine;
        uniform_real_distribution<> uniform;  
        

	private:

        int Niterations;                         // Counter saving the number of nested loops used
        double informationGain;                  // Information gain in moving from prior to posterior PDF
        double logEvidence;                      // Skilling's evidence
        double logEvidenceError;                 // Skilling's error on evidence (based on information gain)
        double logMeanEvidence;                  // Keeton's mean evidence (derived from Skilling equation)
        double logMeanTotalEvidence;             // Keeton's mean live evidence
        double logMeanEvidenceError;             // Keeton's error on evidence (no contribution from live evidence)
        double logMeanTotalEvidenceError;        // Keeton's total error on evidence (contribution from live evidence)
        double logMeanLikelihoodOfLivePoints;    // The logarithm of the mean likelihood value of the remaining set of live points
        double logMaximumLikelihoodOfLivePoints; // The logarithm of the maximun likelihood value of the remaining set of live points
        double computationalTime;
        double constant1;                        // Constant factors used in Keeton's formulas
        double constant2;
        double constant3;
        ArrayXd logLikelihood;                   // log-likelihood values of the actual set of live points
        ArrayXXd nestedSampleOfParameters;       // parameters values (the free parameters of the problem) of the actual set of live points

        void computeKeetonEvidenceError(const bool printFlag, const double logMeanLiveEvidence);
        void printComputationalTime(const double startTime);

}; 

#endif
