// Class for computing final results to be saved in output files.
// Created by Enrico Corsaro @ IvS - 22 February 2013
// e-mail: emncorsaro@gmail.com
// Header file "Results.h"
// Implementations contained in "Results.cpp"


#ifndef RESULTS_H
#define RESULTS_H


#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <string>
#include <iomanip>
#include <cassert>
#include <limits>
#include <Eigen/Core>
#include "Functions.h"
#include "File.h"
#include "NestedSampler.h"


using namespace std;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class Results
{

    public:
 
        Results(NestedSampler &nestedSampler);
        ~Results();
        
        void writeParametersToFile(string fileName, string outputFileExtension = ".txt");
        void writeLogLikelihoodToFile(string fileName);
        void writeLogWeightsToFile(string fileName);
        void writeEvidenceInformationToFile(string fileName);
        void writePosteriorProbabilityToFile(string fileName);
        void writeLogEvidenceToFile(string fileName);
        void writeLogMeanLiveEvidenceToFile(string fileName);
        void writeParametersSummaryToFile(string fileName, const double credibleLevel = 68.27, const bool writeMarginalDistribution = true);
        void writeObjectsIdentificationToFile(){};          // TO DO


    protected:

        double marginalDistributionMode;
        ArrayXd parameterValues;
        ArrayXd marginalDistribution;
        ArrayXd parameterValuesRebinned;
        ArrayXd parameterValuesInterpolated;
        ArrayXd marginalDistributionRebinned;
        ArrayXd marginalDistributionInterpolated;


    private:

        NestedSampler &nestedSampler;
       
        ArrayXd posteriorProbability();
        void writeMarginalDistributionToFile(const int parameterNumber);
        ArrayXd computeCredibleLimits(const double credibleLevel, const double skewness, const int NinterpolationsPerBin = 10);
        ArrayXXd parameterEstimation(const double credibleLevel, const bool writeMarginalDistribution);

};
#endif
