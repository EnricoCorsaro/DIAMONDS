// Derived class for K-means clustering algorithm.
// Created by Joris De Ridder @ IvS - February 2013
// e-mail: emncorsaro@gmail.com
// Header file "KmeansClusterer.h"
// Implementations contained in "KmeansClusterer.cpp"


#ifndef KMEANSCLUSTERER_H
#define KMEANSCLUSTERER_H

#include <ctime>
#include <cfloat>
#include <cmath>
#include <random>
#include <limits>
#include <iostream>
#include "Clusterer.h"


using namespace std;


class KmeansClusterer : public Clusterer
{
    public:
    
        KmeansClusterer(Metric &metric, Projector &featureProjector, bool featureProjectionActivated, 
        unsigned int minNclusters, unsigned int maxNclusters, unsigned int Ntrials, double relTolerance);
        
        ~KmeansClusterer();
    
        virtual int cluster(RefArrayXXd sample, vector<int> &optimalClusterIndices, vector<int> &optimalClusterSizes);
  

    protected:
   

    private:
    
        void chooseInitialClusterCenters(RefArrayXXd sample, RefArrayXXd centers);
        bool updateClusterCentersUntilConverged(RefArrayXXd sample, RefArrayXXd centers, 
                                                RefArrayXd clusterSizes, vector<int> &clusterIndices,
                                                double &sumOfDistancesToClosestCenter, double relTolerance);
        double evaluateBICvalue(RefArrayXXd sample, RefArrayXXd centers, RefArrayXd clusterSizes, 
                                vector<int> &clusterIndices);

        unsigned int minNclusters;
        unsigned int maxNclusters;
        unsigned int Ntrials;
        unsigned int Ndimensions;
        unsigned int Npoints;
        unsigned int Nclusters;
        double relTolerance;
        mt19937 engine;

};



#endif
