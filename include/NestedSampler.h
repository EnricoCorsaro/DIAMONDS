// Class for nested sampling inference
// Enrico Corsaro @ IvS - 24 January 2013
// e-mail: emncorsaro@gmail.com
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
#include "File.h"

using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXi> RefArrayXi;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;

class LivePointsReducer;

class NestedSampler
{
    public:

        NestedSampler(const bool printOnTheScreen, const int initialNlivePoints, const int minNlivePoints, vector<Prior*> ptrPriors, 
                      Likelihood &likelihood, Metric &metric, Clusterer &clusterer); 
        ~NestedSampler();
        
        void run(LivePointsReducer &livePointsReducer, const int NinitialIterationsWithoutClustering = 1000, 
                 const int NiterationsWithSameClustering = 50, const int maxNdrawAttempts = 10000, 
                 const double minRatioOfRemainderToCurrentEvidence = 0.05, const int maxNiterations = 0, 
                 string pathPrefix = "");
        
        virtual bool drawWithConstraint(const RefArrayXXd totalSample, const unsigned int Nclusters, const vector<int> &clusterIndices,
                                        const vector<int> &clusterSizes, RefArrayXd drawnPoint, 
                                        double &logLikelihoodOfDrawnPoint, const int maxNdrawAttempts) = 0;
       

        // Define set and get functions

        unsigned int getNiterations();
        unsigned int getNdimensions();
        int getNlivePoints();
        int getInitialNlivePoints();
        int getMinNlivePoints();
        double getLogCumulatedPriorMass();
        double getLogRemainingPriorMass();
        double getRatioOfRemainderToCurrentEvidence();
        double getLogMaxLikelihoodOfLivePoints();
        double getComputationalTime();
        double getTerminationFactor();
        vector<int> getNlivePointsPerIteration();
        ArrayXXd getNestedSample();
        ArrayXd getLogLikelihood();
        
        void setLogEvidence(double newLogEvidence);
        double getLogEvidence();
        
        void setLogEvidenceError(double newLogEvidenceError);
        double getLogEvidenceError();
        
        void setInformationGain(double newInformationGain);
        double getInformationGain();
        
        void setPosteriorSample(ArrayXXd newPosteriorSample);
        ArrayXXd getPosteriorSample();
        
        void setLogLikelihoodOfPosteriorSample(ArrayXd newLogLikelihoodOfPosteriorSample);
        ArrayXd getLogLikelihoodOfPosteriorSample();
        
        void setLogWeightOfPosteriorSample(ArrayXd newLogWeightOfPosteriorSample);
        ArrayXd getLogWeightOfPosteriorSample();
        
        ArrayXd getLogEvidenceOfPosteriorSample();
        ArrayXd getLogMeanLiveEvidenceOfPosteriorSample();
        
        void setOutputPathPrefix(string newOutputPathPrefix);
        string getOutputPathPrefix();
       
        ofstream outputFile;                        // An output file stream to save configuring parameters also from derived classes 


    protected:

        vector<Prior*> ptrPriors;                   // A vector of pointers to objects of class Prior, containing the priors for each parameter
        Likelihood &likelihood;                     // An object of class Likelihood to contain the likelihood used in the Bayesian inference
        Metric &metric;                             // An object of class Metric for the proper metric to adopt in the computation
        Clusterer &clusterer;                       // An object of class Clusterer to contain the cluster algorithm used in the process
        bool printOnTheScreen;                      // A boolean specifying whether we want current results to be printed on the screen 
        unsigned int Ndimensions;                   // Total number of dimensions of the inference
        int NlivePoints;                            // Total number of live points at a given iteration
        int minNlivePoints;                         // Minimum number of live points allowed
        int reducedNdimensions;                     // Number of effective dimensions of the clustering when a feature projector is activated
        double worstLiveLogLikelihood;              // The worst log likelihood value of the current live sample
        double logCumulatedPriorMass;               // The total (cumulated) prior mass at a given nested iteration
        double logRemainingPriorMass;               // The remaining width in prior mass at a given nested iteration (log X)
        double ratioOfRemainderToCurrentEvidence;   // The current ratio of live to cumulated evidence 
        vector<int> NlivePointsPerIteration;        // A vector that stores the number of live points used at each iteration of the nesting process
        
        mt19937 engine;
        virtual bool verifySamplerStatus() = 0; 
        

	private:

        string outputPathPrefix;                 // The path of the directory where all the results have to be saved
        unsigned int Niterations;                // Counter saving the number of nested loops used
        int updatedNlivePoints;                  // The updated number of live points to be used in the next iteration
        int initialNlivePoints;                  // The initial number of live points
        double informationGain;                  // Skilling's Information gain in moving from prior to posterior PDF
        double logEvidence;                      // Skilling's Evidence
        double logEvidenceError;                 // Skilling's error on Evidence (based on IG)
        double logMaxLikelihoodOfLivePoints;     // The maximum log(Likelihood) of the set of live points
        double logMeanLikelihoodOfLivePoints;    // The logarithm of the mean likelihood value of the current set of live points
        double computationalTime;                // Computational time of the process
        double terminationFactor;                // The final value of the stopping condition for the nested process
        ArrayXXd nestedSample;                   // Parameter values (for all the free parameters of the problem) of the current set of live points
        ArrayXd logLikelihood;                   // log-likelihood values of the current set of live points
                                                 // is removed from the sample.
        ArrayXXd posteriorSample;                // Parameter values (for all the free parameters of the problem) in the final posterior sampling
        ArrayXd logLikelihoodOfPosteriorSample;  // log(Likelihood) values corresponding to the posterior sample 
        ArrayXd logWeightOfPosteriorSample;      // log(Weights) = log(Likelihood) + log(dX) corresponding to the posterior sample
        ArrayXd logEvidenceOfPosteriorSample;    // log(Evidence) cumulated at iteration in the nesting process
        ArrayXd logMeanLiveEvidenceOfPosteriorSample; // log(MeanLiveEvidence) remaining at each iteration in the nesting process

        void removeLivePointsFromSample(const vector<int> &indicesOfLivePointsToRemove, 
                                        vector<int> &clusterIndices, vector<int> &clusterSizes);
        void printComputationalTime(const double startTime);
}; 

#endif
