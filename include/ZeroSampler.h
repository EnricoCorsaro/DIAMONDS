// Derived class for building empty nested sampler objects.
// Created by Enrico Corsaro @ IvS - 9 May 2014
// e-mail: emncorsaro@gmail.com
// Header file "ZeroSampler.h"
// Implementation contained in "ZeroSampler.cpp"

#ifndef ZEROSAMPLER_H
#define ZEROSAMPLER_H

#include <algorithm>
#include <Eigen/Dense>
#include "NestedSampler.h"

using namespace std;

class ZeroSampler : public NestedSampler
{

    public:
       
        ZeroSampler(const bool printOnTheScreen, const int initialNlivePoints, const int minNlivePoints, vector<Prior*> ptrPriors, 
                    Likelihood &likelihood, Metric &metric, Clusterer &clusterer); 
        ~ZeroSampler();

        virtual bool drawWithConstraint(const RefArrayXXd totalSample, const unsigned int Nclusters, const vector<int> &clusterIndices,
                                        const vector<int> &clusterSizes, RefArrayXd drawnPoint, 
                                        double &logLikelihoodOfDrawnPoint, const int maxNdrawAttempts) override; 
        
    protected:
      
    private:

};

#endif
