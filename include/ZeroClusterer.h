// Derived class for zero (null) clustering algorithm.
// This clusterer does not perform any clustering and 
// is intented for constructing empty clusterer objects.
// Created by Enrico Corsaro @ IvS - 9 May 2014
// e-mail: emncorsaro@gmail.com
// Header file "ZeroClusterer.h"
// Implementations contained in "ZeroClusterer.cpp"

#ifndef ZEROCLUSTERER_H
#define ZEROCLUSTERER_H

#include <ctime>
#include <cfloat>
#include <cmath>
#include <random>
#include <limits>
#include <iostream>
#include "Clusterer.h"


using namespace std;


class ZeroClusterer : public Clusterer
{
    public:
    
        ZeroClusterer(Metric &metric, Projector &featureProjector, bool featureProjectionActivated);
        ~ZeroClusterer();
    
        virtual int cluster(RefArrayXXd sample, vector<int> &optimalClusterIndices, vector<int> &optimalClusterSizes) override;


    protected:
        
    private:

};




#endif
