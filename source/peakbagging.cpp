// Main code for peak bagging by means of nested sampling analysis
// Created by Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Source code file "peakbagging.cpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "MathExtra.h"
#include "File.h"
#include "NestedSampler.h"

#include "UniformPrior.h"
#include "NormalLikelihood.h"
#include "MonoLorentzianModel.h"


int main(int argc, char *argv[])
{
    unsigned long Nrows;
    int Ncols;
    ArrayXXd data;
    
    if (argc != 2)
    {
        cerr << "Usage: peakbagging <inputFile> <outputFile>" << endl;
        exit(EXIT_FAILURE);
    }

    ifstream inputFile(argv[1]);
    if (!inputFile.good())
    {
        cerr << "Error opening input file" << endl;
        exit(EXIT_FAILURE);
    }

    File::snifFile(inputFile, Nrows, Ncols);
    data = File::arrayFromFile(inputFile, Nrows, Ncols);
   
    // Creating arrays for each data type
    ArrayXd covariates = data.col(0);
    ArrayXd observations = data.col(1);
    ArrayXd uncertainties = data.col(2);

    ofstream outputFile(argv[2]);
    if (!outputFile.good())
    {
        cerr << "Error opening output file" << endl;
        exit(EXIT_FAILURE);
    }
      
    // Choose fundamental parameters for the nested inference process
    int Nobjects = 100;      // Number of objects per nested iteration (usually 100)
    int Ndimensions = 1;        // Number of free parameters (dimensions) of the problem

    // Define boundaries of the free parameters of the problem (should be done with separate routine)
    ArrayXXd parametersBoundaries(Ndimensions,2);
    parametersBoundaries(0,0) = 0.0;
    parametersBoundaries(0,1) = 20.0;

    // First step - Setting Prior distribution and parameter space
    UniformPrior prior(parametersBoundaries, Nobjects);

    // Second step - Setting up a model for the inference problem
    MonoLorentzianModel model(covariates);

    // Third step - Setting up the likelihood function to be used
    NormalLikelihood likelihood(covariates, observations, uncertainties, model);

    // Fourth step - Starting nested sampling process
    NestedSampler nestedSampler(prior, likelihood);
    nestedSampler.run();

    // Save the results in an output file (should be done with separate routine)
    outputFile << "# Parameter value    logLikelihood" << endl;
    outputFile << setiosflags(ios::fixed) << setprecision(8);
    File::arrayToFile(outputFile, nestedSampler.posteriorSample.row(0), nestedSampler.logLikelihoodOfPosteriorSample);
    outputFile.close();
    
    cerr << " Evidence: logZ = " << nestedSampler.getLogEvidence() << " +/- " << nestedSampler.getLogEvidenceError() << endl;
    cerr << " Information Gain = " << nestedSampler.getInformationGain() << endl;
    
    return EXIT_SUCCESS;
}
