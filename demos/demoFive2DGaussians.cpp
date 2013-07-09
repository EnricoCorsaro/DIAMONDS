//
// Compile with: clang++ -o demoFive2DGaussians demoFive2DGaussians.cpp -L../build/ -I ../include/ -l multinest -stdlib=libc++ -std=c++11
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
#include "LorentzianModel.h"
#include "demoFive2DGaussians.h"



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

    int Ndimensions = 2;        // Number of free parameters (dimensions) of the problem
    vector<Prior*> ptrPriors(1);
    ArrayXd parametersMinima(Ndimensions);
    ArrayXd parametersMaxima(Ndimensions);
    parametersMinima << -0.5, -0.5;         // Centroid x direction, Centroid y direction
    parametersMaxima << +0.5, +0.5;
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;
    
    // Set up a dummy model. This won't be used because we're computing
    // the Likelihood directly, but the Likelihood nevertheless expects a model in 
    // its constructor.
    
    LorentzianModel model(covariates);


    // Set up the likelihood function to be used
    
    Multiple2DGaussiansLikelihood likelihood(observations, model);
    

    // Set up the K-means clusterer using an Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 5;
    int maxNclusters = 5;
    int Ntrials = 200;
    double relTolerance = 0.01;

    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 


    // Start nested sampling process
    
    bool printOnTheScreen = true;               // Print results on the screen
    int Nobjects = 500;                         // TODO
    int maxNdrawAttempts = 500;                 // TODO
    int NiterationsBeforeClustering = 10;       // Number of nesting iterations before executing clustering algorithm again
    double initialEnlargementFactor = 2.5;      // TODO
    double shrinkingRate = 0.3;                 // Exponent for remaining prior mass in ellipsoid enlargement factor
    double terminationFactor = 0.01;            // Termination factor for nesting loop


    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriors, likelihood, myMetric, kmeans, 
                                        Nobjects, initialEnlargementFactor, shrinkingRate);
    nestedSampler.run(terminationFactor, NiterationsBeforeClustering, maxNdrawAttempts);


    // Save the results in output files

    Results results(nestedSampler);
    results.writeParametersToFile("demoFive2DGaussians_Parameter");
    results.writeLogLikelihoodToFile("demoFive2DGaussians_LogLikelihood.txt");
    results.writeEvidenceInformationToFile("demoFive2DGaussians_Evidence.txt");
    results.writePosteriorProbabilityToFile("demoFive2DGaussians_Posterior.txt");
    results.writeParametersSummaryToFile("demoFive2DGaussians_ParametersSummary.txt");
 
    // That's it!

    return EXIT_SUCCESS;
}
