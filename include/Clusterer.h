// Abstract base class for clustering algorithms.
// Created by Joris De Ridder & Enrico Corsaro @ IvS - February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "Clusterer.h"
// Implementations contained in "Clusterer.cpp"


#ifndef CLUSTERER_H
#define CLUSTERER_H

#include <vector>
#include <Eigen/Core>
#include "Metric.h"


using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class Clusterer
{
    public:
    
        Clusterer(Metric &metric);
        ~Clusterer(){};
    
        virtual int cluster(RefArrayXXd sample, vector<int> &optimalClusterIndices, vector<int> &optimalClusterSizes) = 0;


    protected:
        
        Metric &metric;
    

    private:

};


#endif
