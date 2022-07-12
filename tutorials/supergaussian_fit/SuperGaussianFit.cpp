// Main code for Bayesian Inference with DIAMONDS
// Created by Enrico Corsaro @ OACT - February 2017
// e-mail: emncorsaro@gmail.com
// Source code file "SuperGaussianModel.cpp"

// To compile in Mac OS: 
// clang++ -o SuperGaussianFit SuperGaussianFit.cpp -L../../build/ -I ../../include/ -l diamonds -stdlib=libc++ -std=c++11 -Wno-deprecated-register
// To compile in Linux OS:
// g++ -o SuperGaussianFit SuperGaussianFit.cpp -L../../build/ -I../../include/ -ldiamonds -std=c++11 

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include "Functions.h"
#include "File.h"
#include "MultiEllipsoidSampler.h"
#include "KmeansClusterer.h"
#include "EuclideanMetric.h"
#include "Prior.h"
#include "UniformPrior.h"
#include "NormalPrior.h"
#include "SuperGaussianModel.h"
#include "NormalLikelihood.h"
#include "FerozReducer.h"
#include "PowerlawReducer.h"
#include "Results.h"
#include "Ellipsoid.h"
#include "PrincipalComponentProjector.h"

int main(int argc, char *argv[])
{

    // Check number of arguments for main function
    
    if (argc != 3)
    {
        cerr << "Usage: ./SuperGaussianFit <data filename> <prior filename>" << endl;
        exit(EXIT_FAILURE);
    }


    // ---------------------------
    // ----- Read input data -----
    // ---------------------------

    unsigned long Nrows;
    int Ndimensions;              // Number of parameters for which prior distributions are defined
    int Ncols;
    ArrayXXd data;

    string baseInputDirName = "";
    string inputFileName(argv[1]);
    string outputPathPrefix = "SuperGaussianFit_";

    ifstream inputFile;
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nrows, Ncols);
    data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    inputFile.close();


    // Creating arrays for each data type
    
    ArrayXd covariates = data.col(0);
    ArrayXd observations = data.col(1);
    ArrayXd uncertainties = data.col(2);    

    
    // -------------------------------------------------------
    // ----- First step. Set up all prior distributions -----
    // -------------------------------------------------------
    
    // Uniform Prior
    unsigned long Nparameters;

    // ---- Read prior hyper parameters for resolved modes -----
    string inputFileNamePrior(argv[2]);
    File::openInputFile(inputFile, inputFileNamePrior);
    File::sniffFile(inputFile, Nparameters, Ncols);
    ArrayXXd hyperParameters;
    Ndimensions = Nparameters;

    if (Ncols == 1)
    {
        Ncols = 2;
        hyperParameters.conservativeResize(Nparameters, Ncols);
    }

    hyperParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();

    ArrayXd hyperParametersMinima = hyperParameters.col(0);
    ArrayXd hyperParametersMaxima = hyperParameters.col(1);

    int NpriorTypes = 1;                                        // Total number of prior types included in the computation
    vector<Prior*> ptrPriors(NpriorTypes);
    ArrayXd parametersMinima(Ndimensions);
    ArrayXd parametersMaxima(Ndimensions);
    parametersMinima << hyperParametersMinima;      // Minima values for the free parameters (free parameter #1, free parameter #2, ..., etc.)
    parametersMaxima << hyperParametersMaxima;      // Maxima values for the free parameters (same order as minima)
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;

    string fullPathHyperParameters = outputPathPrefix + "hyperParametersUniform.txt";       // Print prior hyper parameters as output
    uniformPrior.writeHyperParametersToFile(fullPathHyperParameters);


    // -------------------------------------------------------------------
    // ---- Second step. Set up the models for the inference problem ----- 
    // -------------------------------------------------------------------
    
    SuperGaussianModel model(covariates);      // Super Gaussian function


    // -----------------------------------------------------------------
    // ----- Third step. Set up the likelihood function to be used -----
    // -----------------------------------------------------------------
    
    NormalLikelihood likelihood(observations, uncertainties, model);
    

    // -------------------------------------------------------------------------------
    // ----- Fourth step. Set up the K-means clusterer using an Euclidean metric -----
    // -------------------------------------------------------------------------------

    inputFileName = "Xmeans_configuringParameters.txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);

    if (Nparameters != 2)
    {
        cerr << "Wrong number of input parameters for X-means algorithm." << endl;
        exit(EXIT_FAILURE);
    }

    ArrayXd configuringParameters;
    configuringParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();
    
    int minNclusters = configuringParameters(0);
    int maxNclusters = configuringParameters(1);
    
    if ((minNclusters <= 0) || (maxNclusters <= 0) || (maxNclusters < minNclusters))
    {
        cerr << "Minimum or maximum number of clusters cannot be <= 0, and " << endl;
        cerr << "minimum number of clusters cannot be larger than maximum number of clusters." << endl;
        exit(EXIT_FAILURE);
    }

    int Ntrials = 10;
    double relTolerance = 0.01;

    bool printNdimensions = false;
    PrincipalComponentProjector projector(printNdimensions);
    bool featureProjectionActivated = false;

    EuclideanMetric myMetric;
    KmeansClusterer kmeans(myMetric, projector, featureProjectionActivated, 
                           minNclusters, maxNclusters, Ntrials, relTolerance); 


    // ---------------------------------------------------------------------
    // ----- Sixth step. Configure and start nested sampling inference -----
    // ---------------------------------------------------------------------
    
    inputFileName = "NSMC_configuringParameters.txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);
    configuringParameters.setZero();
    configuringParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();

    if (Nparameters != 8)
    {
        cerr << "Wrong number of input parameters for NSMC algorithm." << endl;
        exit(EXIT_FAILURE);
    }

    bool printOnTheScreen = true;                       // Print results on the screen
    int initialNobjects = configuringParameters(0);     // Initial number of live points 
    int minNobjects = configuringParameters(1);         // Minimum number of live points 
    int maxNdrawAttempts = configuringParameters(2);    // Maximum number of attempts when trying to draw a new sampling point
    int NinitialIterationsWithoutClustering = configuringParameters(3); // The first N iterations, we assume that there is only 1 cluster
    int NiterationsWithSameClustering = configuringParameters(4);       // Clustering is only happening every N iterations.
    double initialEnlargementFraction = configuringParameters(5);   // Fraction by which each axis in an ellipsoid has to be enlarged.
                                                                    // It can be a number >= 0, where 0 means no enlargement.
    double shrinkingRate = configuringParameters(6);        // Exponent for remaining prior mass in ellipsoid enlargement fraction.
                                                            // It is a number between 0 and 1. The smaller the slower the shrinkage
                                                            // of the ellipsoids.
    double terminationFactor = configuringParameters(7);    // Termination factor for nested sampling process.

    
    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriors, likelihood, myMetric, kmeans, 
                                        initialNobjects, minNobjects, initialEnlargementFraction, shrinkingRate);
    
    double tolerance = 1.e2;
    double exponent = 0.4;
    PowerlawReducer livePointsReducer(nestedSampler, tolerance, exponent, terminationFactor);
 
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
    results.writeLogWeightsToFile("logWeights.txt");
    results.writeEvidenceInformationToFile("evidenceInformation.txt");
    results.writePosteriorProbabilityToFile("posteriorDistribution.txt");
    results.writeLogEvidenceToFile("logEvidence.txt");
    results.writeLogMeanLiveEvidenceToFile("logMeanLiveEvidence.txt");

    double credibleLevel = 68.3;
    bool writeMarginalDistributionToFile = true;
    results.writeParametersSummaryToFile("parameterSummary.txt", credibleLevel, writeMarginalDistributionToFile);

    cout << "Process completed." << endl;
    
    return EXIT_SUCCESS;
}
