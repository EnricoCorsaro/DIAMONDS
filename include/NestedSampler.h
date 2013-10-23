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
#include <limits>
#include <algorithm>
#include <Eigen/Dense>
#include "Functions.h"
#include "Prior.h"
#include "Likelihood.h"
#include "Metric.h"
#include "Clusterer.h"
#include "LivePointsReducer.h"


using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXi> RefArrayXi;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;

class LivePointsReducer;

class NestedSampler
{
    public:

        NestedSampler(const bool printOnTheScreen, const int initialNobjects, const int minNobjects, vector<Prior*> ptrPriors, 
                      Likelihood &likelihood, Metric &metric, Clusterer &clusterer); 
        ~NestedSampler();
        
        void run(LivePointsReducer &livePointsReducer, const double maxRatioOfRemainderToCurrentEvidence = 0.05, const int NinitialIterationsWithoutClustering = 100, 
                 const int NiterationsWithSameClustering = 50, const int maxNdrawAttempts = 5000);
        
        virtual bool drawWithConstraint(const RefArrayXXd totalSample, const int Nclusters, const vector<int> &clusterIndices,
                                        const vector<int> &clusterSizes, RefArrayXd drawnPoint, 
                                        double &logLikelihoodOfDrawnPoint, const int maxNdrawAttempts) = 0;
        
        int getNiterations();
        int getNobjects();
        int getMinNobjects();
        double getLogCumulatedPriorMass();
        double getLogRemainingPriorMass();
        double getLogEvidence();
        double getLogEvidenceError();
        double getInformationGain();
        double getLogMaxLikelihoodOfLivePoints();
        double getComputationalTime();
        double getTerminationFactor();
        ArrayXXd getNestedSample();
        ArrayXd getLogLikelihood();
        ArrayXXd getPosteriorSample();
        ArrayXd getLogLikelihoodOfPosteriorSample();
        ArrayXd getLogWeightOfPosteriorSample();


    protected:

        vector<Prior*> ptrPriors;
        Likelihood &likelihood;
        Metric &metric;
        Clusterer &clusterer;
        bool printOnTheScreen;
        int Ndimensions;
        int Nobjects;                           // Total number of live points at a given iteration
        int minNobjects;                        // Minimum number of live points allowed
        double worstLiveLogLikelihood;          // The worst likelihood value of the current live sample
        double logCumulatedPriorMass;           // The total (cumulated) prior mass at a given nested iteration
        double logRemainingPriorMass;           // The remaining width in prior mass at a given nested iteration (log X)
        mt19937 engine;
        

	private:

        int Niterations;                         // Counter saving the number of nested loops used
        double informationGain;                  // Skilling's Information gain in moving from prior to posterior PDF
        double logEvidence;                      // Skilling's Evidence
        double logEvidenceError;                 // Skilling's error on Evidence (based on IG)
        double logMaxLikelihoodOfLivePoints;     // The maximum log(Likelihood) of the set of live points
        double logMeanLikelihoodOfLivePoints;    // The logarithm of the mean likelihood value of the current set of live points
        double computationalTime;                // Computational time of the process
        double terminationFactor;                // The condition that sets the total number of nested iterations
        ArrayXXd nestedSample;                   // Parameters values (for all the free parameters of the problem) of the current set of live points
        ArrayXd logLikelihood;                   // log-likelihood values of the current set of live points
                                                 // is removed from the sample.
        ArrayXXd posteriorSample;                // Parameter values (for all the free parameters of the problem) in the final posterior sampling
        ArrayXd logLikelihoodOfPosteriorSample;  // log(Likelihood) values corresponding to the posterior sample 
        ArrayXd logWeightOfPosteriorSample;      // log(Weights) = log(Likelihood) + log(dX) corresponding to the posterior sample

        bool updateNobjects(double logMaxEvidenceContributionNew, double maxRatioOfRemainderToCurrentEvidence);
        void printComputationalTime(const double startTime);
}; 

#endif
