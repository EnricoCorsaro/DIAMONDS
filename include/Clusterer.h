// Abstract base class for clustering algorithms.
// Created by Joris De Ridder & Enrico Corsaro @ IvS - February 2013
// Edited by Enrico Corsaro @ OACT - October 2018
// e-mail: emncorsaro@gmail.com
// Header file "Clusterer.h"
// Implementations contained in "Clusterer.cpp"


#ifndef CLUSTERER_H
#define CLUSTERER_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "Metric.h"
#include "Projector.h"


using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class Clusterer
{
    public:
    
        Clusterer(Metric &metric, Projector &featureProjector, bool featureProjectionActivated = false);
        ~Clusterer(){};
    
        virtual int cluster(RefArrayXXd sample, vector<int> &optimalClusterIndices, vector<int> &optimalClusterSizes) = 0;
        unsigned int getReducedNdimensions();


    protected:
        
        Metric &metric;
        Projector &featureProjector;
        bool featureProjectionActivated;

    private:

};


#endif
