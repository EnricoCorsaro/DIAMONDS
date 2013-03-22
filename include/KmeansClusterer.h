
#ifndef KMEANSCLUSTERER_H
#define KMEANSCLUSTERER_H

#include <ctime>
#include <cfloat>
#include <cmath>
#include <random>
#include <iostream>
#include "Clusterer.h"


using namespace std;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;


class KmeansClusterer : public Clusterer
{
    public:
    
        KmeansClusterer(Metric &metric, unsigned int minNclusters, unsigned int maxNclusters, unsigned int Ntrials, double relTolerance);
        ~KmeansClusterer();
    
        int cluster(RefArrayXXd sample, vector<int> &clusterIndices);
    
    protected:
    
    private:
    
        void chooseInitialClusterCenters(RefArrayXXd sample, RefArrayXXd centers, unsigned int Nclusters);
        bool updateClusterCentersUntilConverged(RefArrayXXd sample, RefArrayXXd centers, 
                                                RefArrayXd clusterSizes, vector<int> &clusterIndices,
                                                double &sumOfDistancesToClosestCenter, double relTolerance);
        double evaluateBICvalue(RefArrayXXd sample, RefArrayXXd centers, RefArrayXd clusterSizes, 
                                vector<int> &clusterIndices);

        unsigned int minNclusters;
        unsigned int maxNclusters;
        unsigned int Ntrials;
        mt19937 engine;
        double relTolerance;

};




#endif