// Compile with:
// clang++ -o demoKmeansClusterer demoKmeansClusterer.cpp -L../build/ -I ../include/ -l multinest -stdlib=libc++ -std=c++11
//

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include "File.h"
#include "EuclideanMetric.h"
#include "KmeansClusterer.h"

using namespace std;
using namespace Eigen;


int main()
{
    // Open the input file and read the data (synthetic sampling of a 2D parameter space)
    
    ifstream inputFile;
    File::openInputFile(inputFile, "kmeans_testsample.txt");
    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
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
    vector<int> clusterIndices(Nrows);
    vector<int> clusterSizes;

    optimalNclusters = kmeans.cluster(sample, clusterIndices, clusterSizes, true);
    
    
    // Output the results 
    
    cerr << "Input number of clusters: 5" << endl; 
    cerr << "Optimal number of clusters: " << optimalNclusters << endl;
    // cout << clusterIndices << endl;
    
    int Nclusters = optimalNclusters;
    int Ndimensions = Ncols;
    int Nobjects = Nrows;
    
    // That's it!
 
    return EXIT_SUCCESS;
}
