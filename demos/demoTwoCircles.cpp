//
// Compile with: clang++ -o demoTwoCircles demoTwoCircles.cpp -L../build/ -I ../include/ -l multinest -stdlib=libc++ -std=c++11
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
#include "demoTwoCircles.h"



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

    int Ndimensions = 2;                      // Number of free parameters (dimensions) of the problem
    vector<Prior*> ptrPriors(1);              // One prior, covering both coordinates
    ArrayXd parametersMinima(Ndimensions);
    ArrayXd parametersMaxima(Ndimensions);
    parametersMinima << -7.0, -4.0;         // Centroid x direction, Centroid y direction
    parametersMaxima << +7.0, +4.0;
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;
    
    // Set up a dummy model. This won't be used because we're computing
    // the Likelihood directly, but the Likelihood nevertheless expects a model in 
    // its constructor.
    
    ZeroModel model(covariates);


    // Set up the likelihood function to be used
    
    TwoCirclesLikelihood likelihood(observations, model);
    

    // Set up the K-means clusterer using an Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 1;
    int maxNclusters = 5;
    int Ntrials = 10;
    double relTolerance = 0.01;

    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 


    // Start nested sampling process
    
    bool printOnTheScreen = true;                   // Print results on the screen
    int Nobjects = 400;                             // TODO
    int maxNdrawAttempts = 200;                     // TODO
    int NinitialIterationsWithoutClustering = 100;  // The first N iterations, we assume that there is only 1 cluster
    int NiterationsWithSameClustering = 10;         // Clustering is only happening every X iterations.
    double initialEnlargementFactor = 2.5;          // TODO
    double shrinkingRate = 0.5;                     // Exponent for remaining prior mass in ellipsoid enlargement factor
    double terminationFactor = 0.01;                // Termination factor for nesting loop


    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriors, likelihood, myMetric, kmeans, 
                                        Nobjects, initialEnlargementFactor, shrinkingRate);
    nestedSampler.run(terminationFactor, NinitialIterationsWithoutClustering, NiterationsWithSameClustering, maxNdrawAttempts);


    // Save the results in output files

    Results results(nestedSampler);
    results.writeParametersToFile("demoTwoCircles_Parameter");
    results.writeLogLikelihoodToFile("demoTwoCircles_LogLikelihood.txt");
    results.writeEvidenceInformationToFile("demoTwoCircles_Evidence.txt");
    results.writePosteriorProbabilityToFile("demoTwoCircles_Posterior.txt");
    results.writeParametersSummaryToFile("demoTwoCircles_ParametersSummary.txt");
 
    // That's it!

    return EXIT_SUCCESS;
}
