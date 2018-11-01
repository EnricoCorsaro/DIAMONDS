// Compile with:
// clang++ -o demoEllipsoids7D demoEllipsoids7D.cpp -L../build/ -I ../include/ -l diamonds -stdlib=libc++ -std=c++11 -Wno-deprecated-register
//

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <unordered_set>
#include <Eigen/Core>
#include "File.h"
#include "EuclideanMetric.h"
#include "KmeansClusterer.h"
#include "Ellipsoid.h"
#include "PrincipalComponentProjector.h"

using namespace std;
using namespace Eigen;


int main()
{
    // ------ IDENTIFY CLUSTERS FROM INPUT SAMPLE ------
    // Open the input file and read the data (synthetic sampling of a 5D parameter space)
    
    ifstream inputFile;
    File::openInputFile(inputFile, "pbagging_testsample7D.txt");
    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
    ArrayXXd data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    ArrayXXd sample = data.transpose();
    inputFile.close();


    // Set up the K-means clusterer using a Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 1;
    int maxNclusters = 20;
    int Ntrials = 10;
    double relTolerance = 0.01;

    bool printNdimensions = false;
    PrincipalComponentProjector projector(printNdimensions);
    bool featureProjectionActivated = true;

    KmeansClusterer kmeans(myMetric, projector, featureProjectionActivated, 
                           minNclusters, maxNclusters, Ntrials, relTolerance); 

 
    // Do the clustering, and get for each point the index of the cluster it belongs to

    int optimalNclusters;
    vector<int> clusterIndices(Nrows);
    vector<int> clusterSizes;

    optimalNclusters = kmeans.cluster(sample, clusterIndices, clusterSizes);
    int Nclusters = optimalNclusters; 
   

    // Output the results 
    
    cerr << "Optimal number of clusters: " << optimalNclusters << endl;
    

    // ------ Compute Ellipsoids ------
    
    int Ndimensions = Ncols;
    assert(sample.cols() == clusterIndices.size());
    assert(sample.cols() >= Ndimensions + 1);            // At least Ndimensions + 1 points are required.


    // The enlargement fraction (it is the fraction by which each axis of an ellipsoid is enlarged)

    double enlargementFraction = 0.0;  
    
    
    // Compute "sorted indices" such that clusterIndices[sortedindices[k]] <= clusterIndices[sortedIndices[k+1]]

    vector<int> sortedIndices = Functions::argsort(clusterIndices);


    // beginIndex will take values such that the indices for one particular cluster (# n) will be in 
    // sortedIndex[beginIndex, ..., beginIndex + clusterSize[n] - 1]      

    int beginIndex = 0;


    // Clear whatever was in the ellipsoids collection

    vector<Ellipsoid> ellipsoids;
    ellipsoids.clear();


    // Create an Ellipsoid for each cluster (provided it's large enough)

    for (int i = 0; i < Nclusters; i++)
    {   
        // Skip cluster if number of points is not large enough

        if (clusterSizes[i] < Ndimensions + 1) 
        {
            // Move the beginIndex up to the next cluster

            beginIndex += clusterSizes[i];


            // Continue with the next cluster

            continue;
        }
        else
        {
            // The cluster is indeed large enough to compute an Ellipsoid.

            // Copy those points that belong to the current cluster in a separate Array
            // This is because Ellipsoid needs a contiguous array of points.

            ArrayXXd sampleOfOneCluster(Ndimensions, clusterSizes[i]);

            for (int n = 0; n < clusterSizes[i]; ++n)
            {
                sampleOfOneCluster.col(n) = sample.col(sortedIndices[beginIndex+n]);
            }


            // Move the beginIndex up to the next cluster

            beginIndex += clusterSizes[i];


            // Add ellipsoid at the end of our vector

            ellipsoids.push_back(Ellipsoid(sampleOfOneCluster, enlargementFraction));
        }
    }

    int Nellipsoids = ellipsoids.size();
    cerr << "Nellispids: " << Nellipsoids << endl;
   
    
    // Find which ellipsoids are overlapping and which are not
    
    vector<unordered_set<int>> overlappingEllipsoidsIndices;


    // Remove whatever was in the container before

    overlappingEllipsoidsIndices.clear();
   

    // Make sure that the indices container has the right size

    overlappingEllipsoidsIndices.resize(ellipsoids.size());


    // If Ellipsoid i overlaps with ellipsoid j, than of course ellipsoid j also overlaps with i.
    // The indices are kept in an unordered_set<> which automatically takes care
    // that there are no duplicates.  

    bool ellipsoidMatrixDecompositionIsSuccessful;

    for (int i = 0; i < Nellipsoids-1; ++i)
    {
        for (int j = i+1; j < Nellipsoids; ++j)
        {
            if (ellipsoids[i].overlapsWith(ellipsoids[j], ellipsoidMatrixDecompositionIsSuccessful))
            {
                overlappingEllipsoidsIndices[i].insert(j);
                overlappingEllipsoidsIndices[j].insert(i);
            }
        }
    }
               
    mt19937 engine;
    clock_t clockticks = clock();
    engine.seed(clockticks);
    uniform_real_distribution<> uniform(0.0, 1.0);  
    

    // Get the hyper-volume for each of the ellipsoids and normalize it 
    // to the sum of the hyper-volumes over all the ellipsoids

    vector<double> normalizedHyperVolumes(Nellipsoids);
    
    for (int n=0; n < Nellipsoids; ++n)
    {
        normalizedHyperVolumes[n] = ellipsoids[n].getHyperVolume();
    }

    double sumOfHyperVolumes = accumulate(normalizedHyperVolumes.begin(), normalizedHyperVolumes.end(), 0.0, plus<double>());

    cout << "Nrmalized Hyper-Volumes" << endl;
    for (int n=0; n < Nellipsoids; ++n)
    {
        normalizedHyperVolumes[n] /= sumOfHyperVolumes;
        cout << "Ellipsoid #" << n << "   " << normalizedHyperVolumes[n] << endl;
    }


    // Pick an ellipsoid with a probability according to its normalized hyper-volume
    // First generate a uniform random number between 0 and 1

    double uniformNumber = uniform(engine);


    // Select the ellipsoid that makes the cumulative hyper-volume greater than this random
    // number. Those ellipsoids with a larger hyper-volume will have a greater probability to 
    // be chosen.

    double cumulativeHyperVolume = normalizedHyperVolumes[0];
    int indexOfSelectedEllipsoid = 0;
    
    while (cumulativeHyperVolume < uniformNumber)
    {
        indexOfSelectedEllipsoid++;
        cumulativeHyperVolume += normalizedHyperVolumes[indexOfSelectedEllipsoid];
    }

    cout << "Selected Ellipsoid #: " << indexOfSelectedEllipsoid << endl;
    
    if (overlappingEllipsoidsIndices[indexOfSelectedEllipsoid].empty())
    {
        cout << "No overlaps occurring for the selected ellipsoid" << endl;
    }
    
    // ------ Draw points from the Ellipsoid ------

    int Npoints = 1000;
    ArrayXXd sampleOfDrawnPoints(Npoints,Ndimensions);
    ArrayXd drawnPoint(Ndimensions);

    for (int i=0; i < Npoints; ++i)
    {
        ellipsoids[indexOfSelectedEllipsoid].drawPoint(drawnPoint);
        sampleOfDrawnPoints.row(i) = drawnPoint.transpose();
    }

    ofstream outputFile;
    File::openOutputFile(outputFile,"ellipsoidSample7D.txt");
    File::arrayXXdToFile(outputFile, sampleOfDrawnPoints);

    return EXIT_SUCCESS;
}
