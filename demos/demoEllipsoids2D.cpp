// Compile with:
// clang++ -o demoEllipsoids2D demoEllipsoids2D.cpp -L../build/ -I ../include/ -l diamonds -stdlib=libc++ -std=c++11 -Wno-deprecated-register
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
    // Open the input file and read the data (synthetic sampling of a 2D parameter space)
    
    ifstream inputFile;
    File::openInputFile(inputFile, "kmeans_testsample2D.txt");
    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
    ArrayXXd data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    ArrayXXd sample = data.transpose();
    inputFile.close();


    
    
    // Set up the clusterer using a Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 2;
    int maxNclusters = 10;
    int Ntrials = 20;
    double relTolerance = 0.01;

    bool printNdimensions = false;
    PrincipalComponentProjector projector(printNdimensions);
    bool featureProjectionActivated = false;

    KmeansClusterer clusterer(myMetric, projector, featureProjectionActivated, 
                           minNclusters, maxNclusters, Ntrials, relTolerance); 


    // Do the clustering, and get for each point the index of the cluster it belongs to

    int optimalNclusters;
    vector<int> clusterIndices(Nrows);
    vector<int> clusterSizes;

    optimalNclusters = clusterer.cluster(sample, clusterIndices, clusterSizes);
    int Nclusters = optimalNclusters; 
   

    // Output the results 
    
    cerr << "Input number of clusters: 5" << endl; 
    cerr << "Optimal number of clusters: " << optimalNclusters << endl;
    

    // ------ Compute Ellipsoids ------
    
    int Ndimensions = Ncols;
    assert(sample.cols() == clusterIndices.size());
    assert(sample.cols() >= Ndimensions + 1);            // At least Ndimensions + 1 points are required.


    // The enlargement fraction (it is the fraction by which each axis of an ellipsoid is enlarged)

    double enlargementFraction = 3.0;  
    
    
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

    ArrayXd centerCoordinate(2);
   
    cout << "Nrmalized Hyper-Volumes" << endl;
    for (int n = 0; n < Nellipsoids; ++n)
    {
        normalizedHyperVolumes[n] /= sumOfHyperVolumes;
        centerCoordinate = ellipsoids[n].getCenterCoordinates();
        cerr << "Ellipsoid #" << n << "   " << normalizedHyperVolumes[n] << endl;
        cerr << "Center Coordinates: " << centerCoordinate.transpose() << endl;
        cerr << endl;
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

    cerr << "Selected Ellipsoid #: " << indexOfSelectedEllipsoid << endl;
    centerCoordinate = ellipsoids[indexOfSelectedEllipsoid].getCenterCoordinates();
    cerr << "Center Coordinates: " << centerCoordinate.transpose() << endl;
    cerr << endl;
   

    // ------ Draw points from the Ellipsoid ------

    int Npoints = 10000;    
    ArrayXXd sampleOfDrawnPoints(Npoints,Ndimensions);
    ArrayXd drawnPoint(Ndimensions);

    for (int i=0; i < Npoints; ++i)
    {
        bool newPointIsFound = false;
        
        while (newPointIsFound == false)
        {
            // Draw a new point inside the ellipsoid
            
            ellipsoids[indexOfSelectedEllipsoid].drawPoint(drawnPoint);
            
            
            // Check if the new point is also in other ellipsoids. If the point happens to be 
            // in N overlapping ellipsoids, then accept it only with a probability 1/N. If we
            // wouldn't do this, the overlapping regions in the ellipsoids would be oversampled.

            if (!overlappingEllipsoidsIndices[indexOfSelectedEllipsoid].empty())
            {
                // There are overlaps, so count the number of ellipsoids to which the new
                // point belongs
            
                int NenclosingEllipsoids = 1;

                for (auto index = overlappingEllipsoidsIndices[indexOfSelectedEllipsoid].begin();
                          index != overlappingEllipsoidsIndices[indexOfSelectedEllipsoid].end();
                        ++index)
                {
                    if (ellipsoids[*index].containsPoint(drawnPoint))  
                    {
                        NenclosingEllipsoids = static_cast<int>(DBL_MAX);
                        //NenclosingEllipsoids++;
                    }
                }


                // Only accept the new point with a probability = 1/NenclosingEllipsoids. 
                // If it's not accepted, go immediately back to the beginning of the while loop, 
                // and draw a new point inside the ellipsoid.

                uniformNumber = uniform(engine);
                newPointIsFound = (uniformNumber < 1./NenclosingEllipsoids);
            }
            else
            {
                // There are no ellipsoids overlapping with the selected one, so the point
                // is automatically accepted

                newPointIsFound = true;
            }
        }

        sampleOfDrawnPoints.row(i) = drawnPoint.transpose();
    }

    ofstream outputFile;
    File::openOutputFile(outputFile,"ellipsoidSample2D.txt");
    File::arrayXXdToFile(outputFile, sampleOfDrawnPoints);

    return EXIT_SUCCESS;
}
