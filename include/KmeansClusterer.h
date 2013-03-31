
#ifndef KMEANSCLUSTERER_H
#define KMEANSCLUSTERER_H

#include <ctime>
#include <cfloat>
#include <cmath>
#include <random>
#include <iostream>
#include "Clusterer.h"


using namespace std;


class KmeansClusterer : public Clusterer
{
    public:
    
        KmeansClusterer(Metric &metric, unsigned int minNclusters, unsigned int maxNclusters, unsigned int Ntrials, double relTolerance);
        ~KmeansClusterer();
    
        int cluster(RefArrayXXd sample, RefArrayXi clusterIndices);
    
    protected:
    
    private:
    
        void chooseInitialClusterCenters(RefArrayXXd sample, RefArrayXXd centers, unsigned int Nclusters);
        bool updateClusterCentersUntilConverged(RefArrayXXd sample, RefArrayXXd centers, 
                                                RefArrayXd clusterSizes, RefArrayXi clusterIndices,
                                                double &sumOfDistancesToClosestCenter, double relTolerance);
        double evaluateBICvalue(RefArrayXXd sample, RefArrayXXd centers, RefArrayXd clusterSizes, 
                                RefArrayXi clusterIndices);

        unsigned int minNclusters;
        unsigned int maxNclusters;
        unsigned int Ntrials;
        double relTolerance;
        mt19937 engine;

};




#endif
