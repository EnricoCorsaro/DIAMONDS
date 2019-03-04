//
// Compile with: 
// clang++ -o demoRastriginFunction demoRastriginFunction.cpp -L../build/ -I ../include/ -l diamonds -stdlib=libc++ -std=c++11 -Wno-deprecated-register
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
#include "PowerlawReducer.h"
#include "demoRastriginFunction.h"
#include "PrincipalComponentProjector.h"


int main(int argc, char *argv[])
{
    unsigned long Nrows;
    int Ncols;
    ArrayXXd data;
  
    // Creating dummy arrays for the covariates and the observations.
    // They're not used because we compute our Likelihood directly. 

    ArrayXd covariates;
    ArrayXd observations;


    // -------------------------------------------------------------------
    // ----- First step. Set up the models for the inference problem ----- 
    // ------------------------------------------------------------------- 

    // Set up a dummy model. This won't be used because we're computing
    // the Likelihood directly, but the Likelihood nevertheless expects a model in 
    // its constructor.
    
    ZeroModel model(covariates);


    // -------------------------------------------------------
    // ----- Second step. Set up all prior distributions -----
    // -------------------------------------------------------

    int Ndimensions = 2;        // Number of free parameters (dimensions) of the problem
    vector<Prior*> ptrPriors(1);
    ArrayXd parametersMinima(Ndimensions);
    ArrayXd parametersMaxima(Ndimensions);
    parametersMinima <<  -5.12, -5.12;
    parametersMaxima << 5.12, 5.12;
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;   


    // -----------------------------------------------------------------
    // ----- Third step. Set up the likelihood function to be used -----
    // -----------------------------------------------------------------

    RastriginLikelihood likelihood(observations, model);
    

    // -------------------------------------------------------------------------------
    // ----- Fourth step. Set up the K-means clusterer using an Euclidean metric -----
    // -------------------------------------------------------------------------------

    EuclideanMetric myMetric;
    int minNclusters = 1;
    int maxNclusters = 10;
    int Ntrials = 10;
    double relTolerance = 0.01;

    bool printNdimensions = false;
    PrincipalComponentProjector projector(printNdimensions);
    bool featureProjectionActivated = false;

    KmeansClusterer kmeans(myMetric, projector, featureProjectionActivated, 
                           minNclusters, maxNclusters, Ntrials, relTolerance); 


    // ---------------------------------------------------------------------
    // ----- Sixth step. Configure and start nested sampling inference -----
    // ---------------------------------------------------------------------
    
    bool printOnTheScreen = true;                   // Print results on the screen
    int initialNobjects = 1000;                      // Initial number of active points evolving within the nested sampling process.
    int minNobjects = 400;                          // Minimum number of active points allowed in the nesting process.
    int maxNdrawAttempts = 50000;                    // Maximum number of attempts when trying to draw a new sampling point.
    int NinitialIterationsWithoutClustering = 500;  // The first N iterations, we assume that there is only 1 cluster.
    int NiterationsWithSameClustering = 20;         // Clustering is only happening every X iterations.
    double initialEnlargementFraction = 3.0;        // Fraction by which each axis in an ellipsoid has to be enlarged.
                                                    // It can be a number >= 0, where 0 means no enlargement.
    double shrinkingRate = 0.05;                     // Exponent for remaining prior mass in ellipsoid enlargement fraction.
                                                    // It is a number between 0 and 1. The smaller the slower the shrinkage
                                                    // of the ellipsoids.
    double terminationFactor = 0.05;                // Termination factor for nesting loop.


    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriors, likelihood, myMetric, kmeans, 
                                        initialNobjects, minNobjects, initialEnlargementFraction, shrinkingRate);
    
    double tolerance = 1.e2;
    double exponent = 0.4;
    PowerlawReducer livePointsReducer(nestedSampler, tolerance, exponent, terminationFactor);  

    string outputPathPrefix = "demoRastriginFunction_";
    nestedSampler.run(livePointsReducer, NinitialIterationsWithoutClustering, NiterationsWithSameClustering, 
                      maxNdrawAttempts, terminationFactor, 0, outputPathPrefix);

    nestedSampler.outputFile << "# List of configuring parameters used for the ellipsoidal sampler and X-means" << endl;
    nestedSampler.outputFile << "# Row #1: Minimum Nclusters" << endl;
    nestedSampler.outputFile << "# Row #2: Maximum Nclusters" << endl;
    nestedSampler.outputFile << "# Row #3: Initial Enlargement Fraction" << endl;
    nestedSampler.outputFile << "# Row #4: Shrinking Rate" << endl;
    nestedSampler.outputFile << minNclusters << endl;
    nestedSampler.outputFile << maxNclusters << endl;
    nestedSampler.outputFile << initialEnlargementFraction << endl;
    nestedSampler.outputFile << shrinkingRate << endl;
    nestedSampler.outputFile.close();


    // -------------------------------------------------------
    // ----- Last step. Save the results in output files -----
    // -------------------------------------------------------

    Results results(nestedSampler);
    results.writeParametersToFile("parameter");
    results.writeLogLikelihoodToFile("logLikelihood.txt");
    results.writeEvidenceInformationToFile("evidenceInformation.txt");
    results.writePosteriorProbabilityToFile("posteriorDistribution.txt");

    double credibleLevel = 68.3;
    bool writeMarginalDistributionToFile = true;
    results.writeParametersSummaryToFile("parameterSummary.txt", credibleLevel, writeMarginalDistributionToFile);

 
    // That's it!

    return EXIT_SUCCESS;
}
