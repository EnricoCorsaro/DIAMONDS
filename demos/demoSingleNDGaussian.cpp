//
// Compile with: clang++ -o demoSingleNDGaussian demoSingleNDGaussian.cpp -L../build/ -I ../include/ -l multinest -stdlib=libc++ -std=c++11
// 

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
#include "Results.h"
#include "Ellipsoid.h"
#include "ZeroModel.h"
#include "FerozReducer.h"
#include "demoSingleNDGaussian.h"



int main(int argc, char *argv[])
{
    unsigned long Nrows;
    int Ncols;
    ArrayXXd data;

  
    // Creating dummy arrays for the covariates and the observations.
    // They're not used because we compute our Likelihood directly. 

    ArrayXd covariates;
    ArrayXd observations;


    // Setting Prior distribution and parameter space

    int Ndimensions = 6;        // Number of free parameters (dimensions) of the problem
    vector<Prior*> ptrPriors(1);
    ArrayXd parametersMinima(Ndimensions);
    ArrayXd parametersMaxima(Ndimensions);
    parametersMinima.fill(-20);         
    parametersMaxima.fill(20);
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;
    

    // Set up a dummy model. This won't be used because we're computing
    // the Likelihood directly, but the Likelihood nevertheless expects a model in 
    // its constructor.
    
    ZeroModel model(covariates);


    // Set up the likelihood function to be used
    
    SingleNDGaussianLikelihood likelihood(observations, model, Ndimensions);


    // Set up the K-means clusterer using an Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 1;
    int maxNclusters = 10;
    int Ntrials = 10;
    double relTolerance = 0.01;

    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 


    // Configure nested sampling
    
    bool printOnTheScreen = true;                   // Print results on the screen
    int initialNobjects = 10000;                     // Initial number of active points evolving within the nested sampling process.
    int minNobjects = 500;                          // Minimum number of active points allowed in the nesting process.
    int maxNdrawAttempts = 5000;                    // Maximum number of attempts when trying to draw a new sampling point.
    int NinitialIterationsWithoutClustering = 1000; // The first N iterations, we assume that there is only 1 cluster.
    int NiterationsWithSameClustering = 50;         // Clustering is only happening every X iterations.
    double initialEnlargementFraction = 2.0;        // Fraction by which each axis in an ellipsoid has to be enlarged.
                                                    // It can be a number >= 0, where 0 means no enlargement.
    double shrinkingRate = 0.8;                     // Exponent for remaining prior mass in ellipsoid enlargement fraction.
                                                    // It is a number between 0 and 1. The smaller the slower the shrinkage
                                                    // of the ellipsoids.
    double terminationFactor = 0.05;                // Termination factor for nesting loop.


    // Save configuring parameters into an ASCII file

    ofstream outputFile;
    string fullPath = "demoSingleNDGaussian_configuringParameters.txt";
    File::openOutputFile(outputFile, fullPath);
    File::configuringParametersToFile(outputFile, initialNobjects, minNobjects, inNclusters, maxNclusters, NinitialIterationsWithoutClustering,
                                     NiterationsWithSameClustering, maxNdrawAttempts, initialEnlargementFraction, shrinkingRate, terminationFactor);
    outputFile.close();
   

    // Start the computation

    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriors, likelihood, myMetric, kmeans, 
                                        initialNobjects, minNobjects, initialEnlargementFraction, shrinkingRate);
        
    double toleranceOnEvidence = 0.01;
    FerozReducer ferozReducer(nestedSampler, toleranceOnEvidence);

    nestedSampler.run(terminationFactor, NinitialIterationsWithoutClustering, NiterationsWithSameClustering, maxNdrawAttempts);


    // Save the results in output files

    Results results(nestedSampler);
    results.writeParametersToFile("demoSingleNDGaussian_Parameter");
    results.writeLogLikelihoodToFile("demoSingleNDGaussian_LogLikelihood.txt");
    results.writeEvidenceInformationToFile("demoSingleNDGaussian_Evidence.txt");
    results.writePosteriorProbabilityToFile("demoSingleNDGaussian_Posterior.txt");
    results.writeParametersSummaryToFile("demoSingleNDGaussian_ParametersSummary.txt");


    // That's it!

    return EXIT_SUCCESS;
}
