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
#include "UniformPrior.h"
#include "HyperEllipsoidSampler.h"

using namespace std;
using namespace Eigen;


int main()
{
    // Open the input file and read the data
    
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
    cout << clusterIndices << endl;
    
    // That's it!

    ArrayXd maxima(2);
    ArrayXd minima(2);

    maxima << sample.row(0).maxCoeff(), sample.row(1).maxCoeff();
    minima << sample.row(0).minCoeff(), sample.row(1).minCoeff();

    UniformPrior prior(minima,maxima);
    HyperEllipsoidSampler sampler(prior, myMetric, Nrows, 0.3, 1);

    ArrayXi NpointsPerCluster(1);
    ArrayXXd allClustersCovarianceMatrix(1,1);
    ArrayXd allCentersCoordinates(1); 

    sampler.computeEllipsoids(sample, clusterIndices, allClustersCovarianceMatrix, allCentersCoordinates, NpointsPerCluster);

    cout << NpointsPerCluster << endl;
    return EXIT_SUCCESS;
}
