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
#include "MultiEllipsoidSampler.h"
#include "KmeansClusterer.h"
#include "EuclideanMetric.h"
#include "Prior.h"
#include "UniformPrior.h"
#include "NormalPrior.h"
#include "NormalLikelihood.h"
#include "LorentzianModel.h"
#include "Results.h"
#include "Ellipsoid.h"

int main(int argc, char *argv[])
{
    unsigned long Nrows;
    int Ncols;
    ArrayXXd data;
  

    // Check number of arguments for main function
    
    if (argc != 3)
    {
        cerr << "Usage: peakbagging <input file> <output directory>" << endl;
        exit(EXIT_FAILURE);
    }


    // Read data from input file specified

    ifstream inputFile(argv[1]);
    if (!inputFile.good())
    {
        cerr << "Error opening input file" << endl;
        exit(EXIT_FAILURE);
    }

    File::snifFile(inputFile, Nrows, Ncols);
    data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    inputFile.close();

   
    // Creating arrays for each data type
    
    ArrayXd covariates = data.col(0);
    ArrayXd observations = data.col(1);
    ArrayXd uncertainties = data.col(2);
   

    // Choose fundamental parameters for the nested inference process

    int Nobjects = 100;      // Number of objects per nested iteration (usually 100)
    int Ndimensions = 3;        // Number of free parameters (dimensions) of the problem


    // First step - Setting Prior distribution and parameter space

    ArrayXd parametersMinima(Ndimensions);
    ArrayXd parametersMaxima(Ndimensions);
    parametersMinima <<  0.0, 0.8, 1.0;         // Centroid, Amplitude, Gamma
    parametersMaxima << 20.0, 1.5, 3.0;

    //ArrayXd parametersMean(Ndimensions);
    //ArrayXd parametersSDV(Ndimensions);
    //parametersMean << 9.2,1.2,2.0;
    //parametersSDV << 0.5,0.2,0.2;
    
    vector<Prior*> ptrPriorsVector(1);
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriorsVector[0] = &uniformPrior;
    //NormalPrior normalPrior(parametersMean, parametersSDV);
    //ptrPriorsVector[0] = &normalPrior;


    // Second step - Setting up a model for the inference problem
    
    LorentzianModel model(covariates);
    

    // Third step - Setting up the likelihood function to be used
    
    NormalLikelihood likelihood(covariates, observations, uncertainties, model);
    

    // Fourth step - Set up the K-means clusterer using an Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 2;
    int maxNclusters = 10;
    int Ntrials = 10;
    double relTolerance = 0.01;

    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 


    // Fifth step - Start nested sampling process
    
    int Ndraws = 1;
    int Nobjects = 1000;
    int NiterationsBeforeClustering = 2;        // Number of nesting iterations before executing clustering algorithm again
    double initialEnlargementFactor = 1.2;  
    double alpha = 1.0;                         // Exponent for remaining prior mass in ellipsoid enlargement factor
    double terminationFactor = 0.5;             // Termination factor for nesting loop

    MultiEllipsoidSampler nestedSampler(ptrPriorsVector, likelihood, myMetric, kmeans, Nobjects, initialEnlargementFactor, alpha);
    nestedSampler.run(terminationFactor, NiterationsBeforeClustering);


    // Save the results in output files

    Results results(nestedSampler);
    string outputDirName(argv[2]);
    results.writeParametersToFile(outputDirName + "/parameter");
    results.writeLogLikelihoodToFile(outputDirName + "/likelihood.txt");
    results.writeEvidenceInformationToFile(outputDirName + "/evidence.txt");
    results.writePosteriorProbabilityToFile(outputDirName + "/posterior.txt");
    results.writeParameterEstimationToFile(outputDirName + "/parameters_estimation.txt");
    
    return EXIT_SUCCESS;
}
