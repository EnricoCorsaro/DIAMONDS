// Compile with:
// clang++ -o demoKmeansClusterer demoKmeansClusterer.cpp ../source/*.cpp -I ../include/ -stdlib=libc++ -std=c++0x
//

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include "File.h"
#include "EuclideanMetric.h"
#include "KmeansClusterer.h"
#include "HyperEllipsoidSampler.h"
#include "HyperEllipsoidIntersector.h"
#include "UniformPrior.h"
#include "NormalLikelihood.h"
#include "LorentzianModel.h"

using namespace std;
using namespace Eigen;


int main()
{
    /* ------ BEGINNING OF X-MEANS CLUSTERING DEMO ----- */
    // Open the input file and read the data (synthetic sampling of a 2D parameter space)
    
    ifstream inputFile;
    File::openInputFile(inputFile, "kmeans_testsample.txt");
    unsigned long Nrows;
    int Ncols;

    File::snifFile(inputFile, Nrows, Ncols);
    ArrayXXd data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    ArrayXXd sample = data.transpose();
    inputFile.close();


    // Set up the K-means clusterer using a Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 2;
    int maxNclusters = 10;
    int Ntrials = 10;
    double relTolerance = 0.01;

    KmeansClusterer kmeans(myMetric, minNclusters, maxNclusters, Ntrials, relTolerance); 

 
    // Do the clustering, and get for each point the index of the cluster it belongs to

    int optimalNclusters;
    ArrayXi clusterIndices(Nrows);
    optimalNclusters = kmeans.cluster(sample, clusterIndices);
    
    
    // Output the results 
    
    cerr << "Optimal number of clusters: " << optimalNclusters << endl;
    // cout << clusterIndices << endl;
    
    int Nclusters = optimalNclusters;
    int Ndimensions = Ncols;
    int Nobjects = Nrows;
    
    // That's it!
    /* ------ END OF X-MEANS CLUSTERING DEMO ----- */


    /* ------ BEGINNING OF ELLIPSOIDAL SAMPLING DEMO ----- */
    // Read synthetic data from input file specified

    File::openInputFile(inputFile, "../data/data2.txt");
    File::snifFile(inputFile, Nrows, Ncols);
    data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    inputFile.close();


    // Creating arrays for each data type
    
    ArrayXd covariates = data.col(0);
    ArrayXd observations = data.col(1);
    ArrayXd uncertainties = data.col(2);


    // Define boundaries of the free parameters of the problem (should be done with separate routine)

    ArrayXd parametersMinima(Ndimensions);
    ArrayXd parametersMaxima(Ndimensions);

    parametersMaxima << sample.row(0).maxCoeff(), sample.row(1).maxCoeff();
    parametersMinima << sample.row(0).minCoeff(), sample.row(1).minCoeff();


    // First step - Setting Prior distribution and parameter space

    UniformPrior uniformPrior(parametersMinima, parametersMaxima);


    // Second step - Setting up a model for the inference problem
    
    LorentzianModel model(covariates);
    

    // Third step - Setting up the likelihood function to be used
    
    NormalLikelihood likelihood(covariates, observations, uncertainties, model);
   

    // Testing of EllipsoidSampler class

    int Ndraws = 1000;
    double initialEnlargementFactor = 2.0;
    double alpha = 1.;
    HyperEllipsoidSampler sampler(uniformPrior, likelihood, myMetric, Nobjects, initialEnlargementFactor, alpha);
    ArrayXXd sampleOfParameters(Ndimensions,Ndraws);


    ofstream outputFile;
    File::openOutputFile(outputFile, "drawnsample.txt");
    
    sampler.drawWithConstraint(sample, Nclusters, clusterIndices, 0, sampleOfParameters);
    
    for (int i=0; i < Ndraws; i++)
    {
        File::arrayXXdToFile(outputFile, sampleOfParameters.col(i).transpose());
    }
    outputFile.close();
    
    ArrayXXd allClustersCovarianceMatrix = sampler.getAllClustersCovarianceMatrix();
    ArrayXd allCentersCoordinates = sampler.getAllCentersCoordinates();
    ArrayXi NobjectsPerCluster = sampler.getNobjectsPerCluster();
    ArrayXd allEnlargedEigenvalues = sampler.getAllEnlargedEigenvalues();
    ArrayXd allEigenvalues = sampler.getAllEigenvalues();
    ArrayXXd allEigenvectorsMatrix = sampler.getAllEigenvectorsMatrix();

    cout << "Number of points per cluster: " << endl;
    cout << NobjectsPerCluster << endl;
    cout << "All centers coordinates: " << endl;
    cout << allCentersCoordinates << endl;
    cout << "Matrix of all original covariance matrices: " << endl;
    cout << allClustersCovarianceMatrix << endl;
    cout << "All eigenvalues: " << endl;
    cout << allEigenvalues << endl;
    cout << "All enlarged eigenvalues: " << endl;
    cout << allEnlargedEigenvalues << endl;

    /* ------ END OF ELLIPSOIDAL SAMPLING DEMO ----- */

    return EXIT_SUCCESS;
}
