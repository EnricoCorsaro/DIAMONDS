//
// Compile with: clang++ -o demoTwo2DGaussians demoTwo2DGaussians.cpp -L../build/ -I ../include/ -l multinest -stdlib=libc++ -std=c++11
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
#include "demoTwo2DGaussians.h"



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
    parametersMinima << -5.0, -5.0;         // Centroid x direction, Centroid y direction
    parametersMaxima << +5.0, +5.0;
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;
   

    // Set up a dummy model. This won't be used because we're computing
    // the Likelihood directly, but the Likelihood nevertheless expects a model in 
    // its constructor.
    
    ZeroModel model(covariates);


    // Set up the likelihood function to be used
    
    Multiple2DGaussiansLikelihood likelihood(observations, model);
    

    // Set up the K-means clusterer using an Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 1;
    int maxNclusters = 3;
    int Ntrials = 200;
    double relTolerance = 0.01;

    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 


    // Configure nested sampling

    bool printOnTheScreen = true;                   // Print results on the screen
    int initialNobjects = 200;                      // Initial number of active points evolving within the nested sampling process.
    int minNobjects = 200;                          // Minimum number of active points allowed in the nesting process.
    int maxNdrawAttempts = 1000;                    // Maximum number of attempts when trying to draw a new sampling point.
    int NinitialIterationsWithoutClustering = 100;  // The first N iterations, we assume that there is only 1 cluster.
    int NiterationsWithSameClustering = 10;         // Clustering is only happening every X iterations.
    double initialEnlargementFraction = 2.0;        // Fraction by which each axis in an ellipsoid has to be enlarged.
                                                    // It can be a number >= 0, where 0 means no enlargement.
    double shrinkingRate = 0.6;                     // Exponent for remaining prior mass in ellipsoid enlargement fraction.
                                                    // It is a number between 0 and 1. The smaller the slower the shrinkage
                                                    // of the ellipsoids.
    double terminationFactor = 0.01;                // Termination factor for nesting loop.


    // Save configuring parameters into an ASCII file

    ofstream outputFile;
    string fullPath = "demoTwo2DGaussians_configuringParameters.txt";
    File::openOutputFile(outputFile, fullPath);
    outputFile << "Initial Nojects: " << initialNobjects << endl;
    outputFile << "Minimum Nobjects: " << minNobjects << endl;
    outputFile << "Minimum Nclusters: " << minNclusters << endl;
    outputFile << "Maximum Nclusters: " << maxNclusters << endl;
    outputFile << "NinitialIterationsWithoutClustering: " << NinitialIterationsWithoutClustering << endl;
    outputFile << "NiterationsWithSameClustering: " << NiterationsWithSameClustering << endl;
    outputFile << "maxNdrawAttempts: " << maxNdrawAttempts << endl;
    outputFile << "Initial EnlargementFraction: " << initialEnlargementFraction << endl;
    outputFile << "Shrinking Rate: " << shrinkingRate << endl;
    outputFile << "terminationFactor: " << terminationFactor << endl;
    outputFile.close();


    // Start the computation

    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriors, likelihood, myMetric, kmeans, 
                                        initialNobjects, minNobjects, initialEnlargementFraction, shrinkingRate);
    nestedSampler.run(terminationFactor, NinitialIterationsWithoutClustering, NiterationsWithSameClustering, maxNdrawAttempts);


    // Save the results in output files

    Results results(nestedSampler);
    results.writeParametersToFile("demoTwo2DGaussians_Parameter");
    results.writeLogLikelihoodToFile("demoTwo2DGaussians_LogLikelihood.txt");
    results.writeEvidenceInformationToFile("demoTwo2DGaussians_Evidence.txt");
    results.writePosteriorProbabilityToFile("demoTwo2DGaussians_Posterior.txt");
    results.writeParametersSummaryToFile("demoTwo2DGaussians_ParametersSummary.txt");


    // That's it!

    return EXIT_SUCCESS;
}
