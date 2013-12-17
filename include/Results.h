// Class for computing final results to be saved in output files.
// Created by Enrico Corsaro @ IvS - 22 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
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
        
        
        // Save results to output Egein Arrays

        ArrayXd posteriorProbability();
        ArrayXXd parameterEstimation(const double credibleLevel = 68.27);
        

        // Write to output files
        
        void writeParametersToFile(string pathPrefix, string outputFileExtension = ".txt");
        void writeLogLikelihoodToFile(string fullPath);
        void writeLogWeightsToFile(string fullPath);
        void writeEvidenceInformationToFile(string fullPath);
        void writePosteriorProbabilityToFile(string fullPath);
        void writeParametersSummaryToFile(string fullPath, const double credibleLevel = 68.27);
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

};
#endif
