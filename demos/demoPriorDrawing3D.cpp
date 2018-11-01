// Compile with:
// clang++ -o demoPriorDrawing3D demoPriorDrawing3D.cpp -L../build/ -I ../include/ -l diamonds -stdlib=libc++ -std=c++11 -Wno-deprecated-register
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
#include "UniformPrior.h"
#include "NormalPrior.h"
#include "SuperGaussianPrior.h"
#include "GridUniformPrior.h"
#include "PrincipalComponentProjector.h"

using namespace std;
using namespace Eigen;


int main()
{
    // ------ IDENTIFY CLUSTERS FROM INPUT SAMPLE ------
    // Open the input file and read the data (synthetic sampling of a 2D parameter space)
    
    ifstream inputFile;
    File::openInputFile(inputFile, "onecluster3D.txt");
    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
    ArrayXXd data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    ArrayXXd sample = data.transpose();
    inputFile.close();


    // Set up the K-means clusterer using a Euclidean metric

    EuclideanMetric myMetric;
    int minNclusters = 1;
    int maxNclusters = 1;
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
    
    cerr << "Input number of clusters: 1" << endl; 
    cerr << "Optimal number of clusters: " << optimalNclusters << endl;
    

    // ------ Compute Ellipsoids ------
    
    int Ndimensions = Ncols;
    assert(sample.cols() == clusterIndices.size());
    assert(sample.cols() >= Ndimensions + 1);            // At least Ndimensions + 1 points are required.


    // The enlargement fraction (it is the fraction by which each axis of an ellipsoid is enlarged)

    double enlargementFraction = 3.00;  
    
    
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

    cerr << "Normalized Hyper-Volumes" << endl;
    ArrayXd centerCoordinate(2);
    ArrayXXd covarianceMatrix(2,2);
    
    for (int n = 0; n < Nellipsoids; ++n)
    {
        normalizedHyperVolumes[n] /= sumOfHyperVolumes;
        centerCoordinate = ellipsoids[n].getCenterCoordinates();
        covarianceMatrix = ellipsoids[n].getCovarianceMatrix();
        cerr << "Ellipsoid #" << n << "   " << normalizedHyperVolumes[n] << endl;
        cerr << "Center Coordinates: " << centerCoordinate.transpose() << endl;
        cerr << "Covariance Matrix: " << endl;
        cerr << covarianceMatrix << endl;
   
        MatrixXd T1 = MatrixXd::Identity(Ndimensions+1,Ndimensions+1);
        T1.bottomLeftCorner(1,Ndimensions) = (-1.0) * centerCoordinate.transpose();
        MatrixXd A = MatrixXd::Zero(Ndimensions+1,Ndimensions+1);
        A(Ndimensions,Ndimensions) = -1;
        A.topLeftCorner(Ndimensions,Ndimensions) = covarianceMatrix.matrix().inverse();
        MatrixXd AT = T1*A*T1.transpose();        // Translating to ellipsoid center
     
        //cerr << "Ellipsoidal Matrix: " << endl;
        //cerr << AT << endl;
        //cerr << endl;
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

    
    cerr << "Selected ellipsoid #: " << indexOfSelectedEllipsoid << endl;
    cerr << endl;



    // ------ Set up prior distributions on each coordinate ------
    
    int Npoints = 10000;    
    ArrayXXd sampleOfDrawnPoints(Npoints,Ndimensions);
    ArrayXd drawnPoint(Ndimensions);
   
    /*      MIX PRIOR       UNIFORM-GRID UNIFORM-UNIFORM
    vector<Prior*> ptrPriors(3);
    ArrayXd parametersMinima(1);
    ArrayXd parametersMaxima(1);
    parametersMinima <<  0.0;
    parametersMaxima << 4.0;
    UniformPrior uniformPrior1(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior1;  

    ArrayXd parametersStartingCoordinate(1);
    ArrayXd parametersNgridPoints(1);
    ArrayXd parametersSeparation(1);
    ArrayXd parametersTolerance(1);
    parametersStartingCoordinate << 0.0;
    parametersNgridPoints << 6;
    parametersSeparation << 0.5;
    parametersTolerance << 0.1;
    GridUniformPrior gridUniformPrior(parametersStartingCoordinate, parametersNgridPoints, parametersSeparation, parametersTolerance);
    ptrPriors[0] = &gridUniformPrior;  

    parametersMinima <<  0.0;
    parametersMaxima << 4.0;
    UniformPrior uniformPrior2(parametersMinima, parametersMaxima);
    ptrPriors[2] = &uniformPrior2;  
    */

    /*      MIX PRIOR       NORMAL-UNIFORM-NORMAL
    vector<Prior*> ptrPriors(3);
    ArrayXd parametersMean(1);
    ArrayXd parametersSDV(1);
    parametersMean <<  2.0;
    parametersSDV << 0.4;
    NormalPrior normalPrior1(parametersMean, parametersSDV);
    ptrPriors[0] = &normalPrior1;  
    
    ArrayXd parametersMinima(1);
    ArrayXd parametersMaxima(1);
    parametersMinima <<  0.0;
    parametersMaxima << 4.0;
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[1] = &uniformPrior;  

    parametersMean <<  2.0;
    parametersSDV << 0.4;
    NormalPrior normalPrior2(parametersMean, parametersSDV);
    ptrPriors[2] = &normalPrior2;  
    */
    
    /*      UNIFORM PRIOR       */
    vector<Prior*> ptrPriors(1);
    ArrayXd parametersMinima(Ndimensions);
    ArrayXd parametersMaxima(Ndimensions);
    parametersMinima <<  0.0, 0.0, 0.0;
    parametersMaxima << 4.0, 4.0, 4.0;
    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;  


    /*      GAUSSIAN PRIOR
    vector<Prior*> ptrPriors(1);
    ArrayXd parametersMean(Ndimensions);
    ArrayXd parametersSDV(Ndimensions);
    parametersMean <<  2.0, 2.0, 2.0;
    parametersSDV << 0.2, 0.4, 0.2;
    NormalPrior normalPrior(parametersMean, parametersSDV);
    ptrPriors[0] = &normalPrior;
    */  
    
    /*      SUPER GAUSSIAN PRIOR
    vector<Prior*> ptrPriors(1);
    ArrayXd parametersMean(Ndimensions);
    ArrayXd parametersSDV(Ndimensions);
    ArrayXd parametersWOP(Ndimensions);
    parametersMean <<  2.0, 2.0, 2.0;
    parametersSDV << 0.1, 0.2, 0.3;
    parametersWOP << 0.4, 0.4, 0.4;
    SuperGaussianPrior superGaussianPrior(parametersMean, parametersSDV, parametersWOP);
    ptrPriors[0] = &superGaussianPrior;  
    */


    // ------ Draw points from the Ellipsoid ------

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
                        //NenclosingEllipsoids = static_cast<int>(DBL_MAX);       // No drawing from overlapping regions!
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


            // The point should not only be drawn inside the ellipsoid, it should also be drawn
            // from the prior. Therefore, accept the point only with the probability given by the
            // prior, so that the regions inside the ellipsoid with a higher prior density will 
            // be sampled more than the regions with a lower prior density. 

            
            // Since different coordinates of our new point may have different priors, 
            // we need to check this for all the priors.
        
            int beginIndex = 0;

            for (int priorIndex = 0; priorIndex < ptrPriors.size(); ++priorIndex)
            {
                // Figure out the number of parameters (=coordinates) that the current prior covers.

                const int NdimensionsOfPrior = ptrPriors[priorIndex]->getNdimensions();


                // Define a subset of the new point, consisting of those coordinates covered by the 
                // same (current) prior distribution.

                ArrayXd subsetOfNewPoint = drawnPoint.segment(beginIndex, NdimensionsOfPrior);


                // Check if the new point is accepted according to the corresponding prior distribution.
             
                newPointIsFound = ptrPriors[priorIndex]->drawnPointIsAccepted(subsetOfNewPoint);
                
                if (!newPointIsFound)
                break;


                // Move the beginIndex on to the next set of coordinates covered by the prior.

                beginIndex += NdimensionsOfPrior;
            }

        }

        sampleOfDrawnPoints.row(i) = drawnPoint.transpose();
    }

    ofstream outputFile;
    File::openOutputFile(outputFile,"priorDrawing3D.txt");
    File::arrayXXdToFile(outputFile, sampleOfDrawnPoints);

    return EXIT_SUCCESS;
}
