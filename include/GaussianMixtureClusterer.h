// Derived class for Gaussian Mixture Model clustering algorithm based on Expectation Maximization.
// Created by Enrico Corsaro @ INAF-OACT - September 2018
// e-mail: emncorsaro@gmail.com
// Header file "GaussianMixtureClusterer.h"
// Implementations contained in "GaussianMixtureClusterer.cpp"


#ifndef GAUSSIANMIXTURECLUSTERER_H
#define GAUSSIANMIXTURECLUSTERER_H

#include <ctime>
#include <cfloat>
#include <cmath>
#include <random>
#include <limits>
#include <iostream>
#include "Clusterer.h"
#include "Functions.h"


using namespace std;


class GaussianMixtureClusterer : public Clusterer
{
    public:
    
        GaussianMixtureClusterer(Metric &metric, Projector &featureProjector, bool featureProjectionActivated,
        unsigned int minNclusters, unsigned int maxNclusters, unsigned int Ntrials, double relTolerance);
        ~GaussianMixtureClusterer();
    
        virtual int cluster(RefArrayXXd sample, vector<int> &optimalClusterIndices, vector<int> &optimalClusterSizes);
        ArrayXXd getCenters();
        ArrayXXd getCovarianceMatrices();
        ArrayXXd getInverseOfCovarianceMatrices();
        ArrayXd getDeterminantOfCovarianceMatrices();


    protected:
  
        ArrayXXd centers;
        ArrayXXd covarianceMatrices;
        ArrayXXd inverseOfCovarianceMatrices;
        ArrayXd determinantOfCovarianceMatrices;
        ArrayXXd differenceFromCenters;
        ArrayXXd multivariateGaussians;
        ArrayXd amplitudes;
        ArrayXd modelProbability;
        ArrayXXd assignmentProbabilities;
        ArrayXd responsibilities;

        double totalLogOfModelProbability;
        double minTotalLogOfModelProbability;

    private:
    
        void chooseInitialClusterCenters(RefArrayXXd sample);
        void chooseInitialClusterCovarianceMatrices(RefArrayXXd sample);
        void computeGaussianMixtureModel(RefArrayXXd sample);
        
        bool updateClustersUntilConverged(RefArrayXXd sample);
        bool searchForEmptyClusters();
        void obtainClusterMembership(vector<int> &optimalClusterIndices, vector<int> &optimalClusterSizes,
                                        RefArrayXXd optimalAssignmentProbabilities);

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
