// Derived class for sampling from an ellipsoid build
// from a sample of points.
// Created by Enrico Corsaro @ IvS - 28 March 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "HyperEllipsoidSampler.h"
// Implementation contained in "HyperEllipsoidSampler.cpp"

#ifndef HYPERELLIPSOIDSAMPLER_H
#define HYPERELLIPSOIDSAMPLER_H

#include <Eigen/Dense>
#include "HyperQuadricSampler.h"
#include "HyperEllipsoidIntersector.h"

using namespace std;

class HyperEllipsoidSampler : public HyperQuadricSampler
{

    public:
       
        HyperEllipsoidSampler(Prior &prior, Likelihood &likelihood, Metric &metric, const int Nobjects, const double initialEnlargementFactor, const double alpha);
        ~HyperEllipsoidSampler();
        
        virtual void drawWithConstraint(const RefArrayXXd totalSampleOfParameters, const int Nclusters, const RefArrayXi clusterIndices,
                                        const double logWidthInPriorMass, RefArrayXXd nestedSampleOfParameters);
        void computeEllipsoids(const RefArrayXXd totalSampleOfParameters, const int Nclusters, 
                               const RefArrayXi clusterIndices, const double logWidthInPrioMass);
        ArrayXXd getAllClustersCovarianceMatrix();
        ArrayXd getAllCentersCoordinates();
        ArrayXi getNobjectsPerCluster();
        ArrayXXd getAllEigenvectorsMatrix();
        ArrayXd getAllEigenvalues();
        ArrayXd getAllEnlargedEigenvalues();
        ArrayXXd getAllEnlargedCovarianceMatrix();
        ArrayXd getHyperVolumes();

    
    protected:
      
        void ellipsoidEnlarger(RefArrayXd eigenvalues, const double logWidthInPriorMass, const int NobjectsInCluster);
        void drawFromHyperSphere(const RefArrayXd eigenvalues, const RefArrayXXd eigenvectorsMatrix, 
                                 const RefArrayXd centerCoordinates, RefArrayXd drawnParameters);
        void computeAllEnlargedCovarianceMatrix();


    private:

        mt19937 engine;
        uniform_real_distribution<> uniform;
        normal_distribution<> normal;
        double initialEnlargementFactor;
        double alpha;
        ArrayXXd allClustersCovarianceMatrix;
        ArrayXXd allEnlargedCovarianceMatrix;
        ArrayXd allCentersCoordinates;
        ArrayXi NobjectsPerCluster;
        ArrayXXd allEigenvectorsMatrix;
        ArrayXd allEigenvalues;
        ArrayXd allEnlargedEigenvalues;
        ArrayXd hyperVolumes;
};

#endif
