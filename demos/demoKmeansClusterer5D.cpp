// Compile with:
// clang++ -o demoKmeansClusterer5D demoKmeansClusterer5D.cpp -L../build/ -I ../include/ -l diamonds -stdlib=libc++ -std=c++11 -Wno-deprecated-register
//

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include "File.h"
#include "EuclideanMetric.h"
#include "KmeansClusterer.h"
#include "PrincipalComponentProjector.h"

using namespace std;
using namespace Eigen;


int main()
{
    // Open the input file and read the data (synthetic sampling of a 5D parameter space)
    
    ifstream inputFile;
    File::openInputFile(inputFile, "kmeans_testsample5D.txt");
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

    bool printNdimensions = true;
    PrincipalComponentProjector projector(printNdimensions);
    bool featureProjectionActivated = false;

    KmeansClusterer kmeans(myMetric, projector, featureProjectionActivated, 
                           minNclusters, maxNclusters, Ntrials, relTolerance); 

 
    // Do the clustering, and get for each point the index of the cluster it belongs to

    int optimalNclusters;
    vector<int> clusterIndices(Nrows);
    vector<int> clusterSizes;

    optimalNclusters = kmeans.cluster(sample, clusterIndices, clusterSizes);
    
    
    // Output the results 
    
    cerr << "Input number of clusters: 4" << endl; 
    cerr << "Optimal number of clusters: " << optimalNclusters << endl;
    
    ArrayXXd finalSample(Nrows,Ncols+1);
    finalSample.block(0,0,Nrows,Ncols) = data;
    
    for (int n = 0; n < Nrows; ++n)
    {
        finalSample(n,Ncols) = clusterIndices[n];
    }

    ofstream outputFile;
    File::openOutputFile(outputFile, "clusterMembershipFromKmeans5D.txt");
    outputFile << scientific << setprecision(4);
    File::arrayXXdToFile(outputFile, finalSample);
    outputFile.close();
    
    // That's it!
 
    return EXIT_SUCCESS;
}
