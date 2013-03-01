// Class for computing final results to be saved in output files.
// Created by Enrico Corsaro @ IvS - 22 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "Results.h"
// Implementations contained in "Results.cpp"


#ifndef RESULTS_H
#define RESULTS_H


#include <cstdlib>
#include <string>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <cassert>
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
 
        Results(NestedSampler &nestedSampler, const char *inputDirectory, const char *inputFilename, const char *outputDirectory);
        ~Results();

        void printParameters();
        void printLogLikelihood();
        void printEvidence();
        void printPosterior();
        void printInference(const double credibleLevel = 68.27);
        ArrayXd getPosteriorDistribution();
        ArrayXXd getInferenceResults();


    private:

        const char *outputDirectory;
        string dataFilename;
        NestedSampler &nestedSampler;
        ArrayXd posteriorDistribution;
        ArrayXXd inferenceResults;

}; // END class Results
#endif
