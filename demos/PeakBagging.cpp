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
#include "ExponentialLikelihood.h"
#include "LorentzianModel.h"
#include "RegularPatternModel.h"
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
   

    // First step - Setting Prior distribution and parameter space

    int Norders = 4;            // Number of radial orders expected in the PSD
    int Ndimensions = 4;        // Number of free parameters (dimensions) of the problem
    vector<Prior*> ptrPriorsVector(1);
    
    ///*
    ArrayXd parametersMinima(Ndimensions);
    ArrayXd parametersMaxima(Ndimensions);
    parametersMinima <<  223.0, 14.0, 1.5, -1.0;         // nuMax, DeltaNu, deltaNu02, deltaNu01
    parametersMaxima << 230.0, 20.0, 4.5, 1.0;
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriorsVector[0] = &uniformPrior;
    //*/

    /*
    ArrayXd parametersMean(Ndimensions);
    ArrayXd parametersSDV(Ndimensions);
    parametersMean << 12.0,1.4,1.5;
    parametersSDV << 2.0,0.5,0.5;
    NormalPrior normalPrior(parametersMean, parametersSDV);
    ptrPriorsVector[0] = &normalPrior;
    */ 


    // Second step - Set up a model for the inference problem
    
    RegularPatternModel model(covariates, Norders);
    

    // Third step - Set up the likelihood function to be used
    
    ExponentialLikelihood likelihood(observations, model);
    

    // Fourth step - Set up the K-means clusterer using an Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 1;
    int maxNclusters = 6;
    int Ntrials = 10;
    double relTolerance = 0.01;

    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 


    // Fifth step - Start nested sampling process
    
    bool printOnTheScreen = true;               // Print results on the screen 
    int NloopMaximum = 5000;                    // Maximum number of attempts when drawing a new active point subject to likelihood constraint
    int Nobjects = 300;                        // Number of active points used in the nesting process
    int NiterationsBeforeClustering = 50;        // Number of nesting iterations before executing clustering algorithm again
    double initialEnlargementFactor = 3.0;      // The initial value for the enlargement factor of the ellipsoids
    double alpha = 0.6;                         // Exponent for shrinkage rate of the ellipsoid enlargement factor
    double terminationFactor = 0.05;            // Termination factor for nesting process

    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriorsVector, 
                                        likelihood, myMetric, kmeans, 
                                        Nobjects, initialEnlargementFactor, alpha);
    nestedSampler.run(terminationFactor, NiterationsBeforeClustering, NloopMaximum);


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
