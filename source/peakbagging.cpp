// Main code for peak bagging by means of nested sampling analysis
// Created by Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Source code file "peakbagging.cpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Functions.h"
#include "File.h"
#include "NestedSampler.h"

#include "UniformPrior.h"
#include "NormalLikelihood.h"
#include "LorentzianModel.h"
#include "Results.h"


int main(int argc, char *argv[])
{
    unsigned long Nrows;
    int Ncols;
    ArrayXXd data;
  

    // Check number of arguments for main function
    
    if (argc != 4)
    {
        cerr << "Usage: peakbagging <input directory> <input data filename> <output directory>" << endl;
        exit(EXIT_FAILURE);
    }


    // Read data from input file specified

    const char *filename = 0;
    string inputFilename;
    string inputDirectory;

    inputFilename = argv[1];
    inputDirectory = argv[2];
    filename = (inputFilename + inputDirectory).data();

    ifstream inputFile(filename);
    if (!inputFile.good())
    {
        cerr << "Error opening input file" << endl;
        exit(EXIT_FAILURE);
    }

    File::snifFile(inputFile, Nrows, Ncols);
    inputFile.clear();
    inputFile.seekg(ios::beg);
    data = File::arrayFromFile(inputFile, Nrows, Ncols);
    inputFile.close();

   
    // Creating arrays for each data type
    
    ArrayXd covariates = data.col(0);
    ArrayXd observations = data.col(1);
    ArrayXd uncertainties = data.col(2);
   

    // Choose fundamental parameters for the nested inference process

    int Nobjects = 100;      // Number of objects per nested iteration (usually 100)
    int Ndimensions = 3;        // Number of free parameters (dimensions) of the problem


    // Define boundaries of the free parameters of the problem (should be done with separate routine)

    ArrayXd parametersMinima(Ndimensions);
    ArrayXd parametersMaxima(Ndimensions);
    parametersMinima(0) = 0.0;
    parametersMaxima(0) = 20.0;
    parametersMinima(1) = 0.8;
    parametersMaxima(1) = 1.5;
    parametersMinima(2) = 1.0;
    parametersMaxima(2) = 3.0;


    // First step - Setting Prior distribution and parameter space

    UniformPrior prior(parametersMinima, parametersMaxima);


    // Second step - Setting up a model for the inference problem
    
    LorentzianModel model(covariates);
    

    // Third step - Setting up the likelihood function to be used
    
    NormalLikelihood likelihood(covariates, observations, uncertainties, model);
    

    // Fourth step - Starting nested sampling process
    
    NestedSampler nestedSampler(prior, likelihood);
    nestedSampler.run(Nobjects);


    // Save the results in an output file

    Results results(nestedSampler, argv[1], argv[2], argv[3]);
    results.printParameters();
    results.printLogLikelihood();
    results.printEvidence();
    results.printPosterior();
    results.printInference();
    
    return EXIT_SUCCESS;
}
