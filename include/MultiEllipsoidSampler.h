// Derived class for sampling from an ellipsoid build
// from a sample of points according to the
// nesting algortihm. The class is based on the
// sampling technique of the MultiNest code presented 
// by Feroz F. et al. (2008, 2009).
// Created by Enrico Corsaro @ IvS - 28 March 2013
// e-mail: emncorsaro@gmail.com
// Header file "MultiEllipsoidSampler.h"
// Implementation contained in "MultiEllipsoidSampler.cpp"

#ifndef MULTIELLIPSOIDSAMPLER_H
#define MULTIELLIPSOIDSAMPLER_H

#include <random>
#include <unordered_set>
#include <algorithm>
#include <Eigen/Dense>
#include "NestedSampler.h"
#include "Ellipsoid.h"

using namespace std;

class MultiEllipsoidSampler : public NestedSampler
{

    public:
       
        MultiEllipsoidSampler(const bool printOnTheScreen, vector<Prior*> ptrPriors, 
                              Likelihood &likelihood, Metric &metric, Clusterer &clusterer, 
                              const int initialNlivePoints, const int minNlivePoints, 
                              const double initialEnlargementFraction, const double shrinkingRate);
        ~MultiEllipsoidSampler();

        virtual bool drawWithConstraint(const RefArrayXXd totalSample, const unsigned int Nclusters, const vector<int> &clusterIndices,
                                        const vector<int> &clusterSizes, RefArrayXd drawnPoint, 
                                        double &logLikelihoodOfDrawnPoint, const int maxNdrawAttempts) override; 
        
        virtual bool verifySamplerStatus() override;

        vector<Ellipsoid> getEllipsoids();
        double getInitialEnlargementFraction();
        double getShrinkingRate();


    protected:
      
        bool ellipsoidMatrixDecompositionIsSuccessful;  // A boolean specifying whether an error occurred in the 
                                                        // eigenvalues decomposition of the ellipsoid matrix
        
        void computeEllipsoids(RefArrayXXd const totalSample, const unsigned int Nclusters, 
                               const vector<int> &clusterIndices, const vector<int> &clusterSizes);
        void findOverlappingEllipsoids(vector<unordered_set<int>> &overlappingEllipsoidsIndices);
        double updateEnlargementFraction(const int clusterSize);


    private:

        vector<Ellipsoid> ellipsoids;
        int Nellipsoids;                        // Total number of ellipsoids computed
        double initialEnlargementFraction;      // Initial fraction for enlargement of ellipsoids
        double shrinkingRate;                   // Prior volume shrinkage rate (between 0 and 1)
        
        uniform_real_distribution<> uniform; 

};

#endif
