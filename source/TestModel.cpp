// Main code for peak bagging by means of nested sampling analysis
// Created by Enrico Corsaro @ IvS - 28 May 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Source code file "TestModel.cpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Functions.h"
#include "File.h"
#include "MultiEllipsoidSampler.h"
#include "KmeansClusterer.h"
#include "EuclideanMetric.h"
#include "TestLikelihood1.h"
#include "TestLikelihood2.h"
#include "Prior.h"
#include "UniformPrior.h"
#include "NormalPrior.h"
#include "Results.h"
#include "Ellipsoid.h"
#include "LorentzianModel.h"

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



    // Setting Prior distribution and parameter space

    int Ndimensions = 2;        // Number of free parameters (dimensions) of the problem
    vector<Prior*> ptrPriorsVector(1);
    
    ///*
    ArrayXd parametersMinima(Ndimensions);
    ArrayXd parametersMaxima(Ndimensions);
    parametersMinima <<  0.0, 0.0;         // Parameter 1, Parameter2
    parametersMaxima << 10.0*Functions::PI, 10.0*Functions::PI;
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


    // Set up a model for the inference problem
    
    LorentzianModel model(covariates);


    // Set up the likelihood function to be used
    
    //TestLikelihood1 likelihood(observations, uncertainties, model);
    TestLikelihood2 likelihood(observations, uncertainties, model);
    

    // Set up the K-means clusterer using an Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 1;
    int maxNclusters = 20;
    int Ntrials = 10;
    double relTolerance = 0.01;

    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 


    // Start nested sampling process
    
    int Nobjects = 2000;
    int NiterationsBeforeClustering = 50;        // Number of nesting iterations before executing clustering algorithm again
    double initialEnlargementFactor = 3.0;  
    double alpha = 0.7;                         // Exponent for remaining prior mass in ellipsoid enlargement factor
    double terminationFactor = 0.05;             // Termination factor for nesting loop

    MultiEllipsoidSampler nestedSampler(ptrPriorsVector, likelihood, myMetric, kmeans, Nobjects, initialEnlargementFactor, alpha);
    nestedSampler.run(terminationFactor, NiterationsBeforeClustering);


    // Save the results in output files

    Results results(nestedSampler);
    string outputDirName(argv[2]);
    results.writeParametersToFile(outputDirName + "/test2_Parameter");
    results.writeLogLikelihoodToFile(outputDirName + "/test2_LogLikelihood.txt");
    results.writeEvidenceInformationToFile(outputDirName + "/test2_Evidence.txt");
    results.writePosteriorProbabilityToFile(outputDirName + "/test2_Posterior.txt");
    
    return EXIT_SUCCESS;
}
