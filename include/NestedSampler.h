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

class NestedSampler
{
    public:

        NestedSampler(const bool printOnTheScreen, const int initialNobjects, const int minNobjects, vector<Prior*> ptrPriors, 
                      Likelihood &likelihood, Metric &metric, Clusterer &clusterer, LivePointsReducer &livePointsReducer); 
        ~NestedSampler();
        
        void run(const double maxRatioOfRemainderToCurrentEvidence = 0.05, const int NinitialIterationsWithoutClustering = 100, 
                 const int NiterationsWithSameClustering = 50, const int maxNdrawAttempts = 5000);

        virtual bool drawWithConstraint(const RefArrayXXd totalSample, const int Nclusters, const vector<int> &clusterIndices,
                                        const vector<int> &clusterSizes, RefArrayXd drawnPoint, 
                                        double &logLikelihoodOfDrawnPoint, const int maxNdrawAttempts) = 0;
        
        int getNiterations();
        double getLogEvidence();
        double getLogEvidenceError();
        double getInformationGain();
        double getComputationalTime();
        ArrayXXd getPosteriorSample();
        ArrayXd getLogLikelihoodOfPosteriorSample();
        ArrayXd getLogWeightOfPosteriorSample();


    protected:

        vector<Prior*> ptrPriors;
        Likelihood &likelihood;
        Metric &metric;
        Clusterer &clusterer;
        LivePointsReducer &livePointsReducer;
        bool printOnTheScreen;
        int Ndimensions;
        int Nobjects;                           // Total number of objects at a given iteration
        double worstLiveLogLikelihood;          // The worst likelihood value of the current live sample
        double logCumulatedPriorMass;           // The total (cumulated) prior mass at a given nested iteration
        double logRemainingPriorMass;           // The remaining width in prior mass at a given nested iteration (log X_k)
        mt19937 engine;
        

	private:

        int minNobjects;                         // Minimum number of live points allowed
        int Niterations;                         // Counter saving the number of nested loops used
        double informationGain;                  // Information gain in moving from prior to posterior PDF
        double logEvidence;                      // Skilling's evidence
        double logEvidenceError;                 // Skilling's error on evidence (based on information gain)
        double logMeanLikelihoodOfLivePoints;    // The logarithm of the mean likelihood value of the actual set of live points
        double logMaxLikelihoodOfLivePoints;     // The logarithm of the maxumum likelihood value of the actual set of live points
        double logMaxEvidenceContribution;       // The logarithm of the maximum evidence contribution at a given iteration of the nesting process
        double computationalTime;                // Computational time of the process
        ArrayXd logLikelihood;                   // log-likelihood values of the actual set of live points
        ArrayXXd nestedSample;                   // parameters values (the free parameters of the problem) of the actual set of live points
        ArrayXXd posteriorSample;                // Parameter values (for all the free parameters of the problem) in the final posterior sampling
        ArrayXd logLikelihoodOfPosteriorSample;  // log(Likelihood) values corresponding to the posterior sample 
        ArrayXd logWeightOfPosteriorSample;      // log(Weights) = log(Likelihood) + log(dX) corresponding to the posterior sample

        bool updateNobjects(double logMaxEvidenceContributionNew, double maxRatioOfRemainderToCurrentEvidence);
        void printComputationalTime(const double startTime);
}; 

#endif
