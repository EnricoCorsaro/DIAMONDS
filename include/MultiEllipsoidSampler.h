// Derived class for sampling from an ellipsoid build
// from a sample of points according to the
// nesting algortihm. The class is based on the
// sampling technique of the MultiNest code presented 
// by Feroz F. et al. (2008, 2009).
// Created by Enrico Corsaro @ IvS - 28 March 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "MultiEllipsoidSampler.h"
// Implementation contained in "MultiEllipsoidSampler.cpp"

#ifndef MULTIELLIPSOIDSAMPLER_H
#define MULTIELLIPSOIDSAMPLER_H

#include <Eigen/Dense>
#include "NestedSampler.h"
#include "Ellipsoid.h"

using namespace std;

class MultiEllipsoidSampler : public NestedSampler
{

    public:
       
        MultiEllipsoidSampler(Prior &prior, Likelihood &likelihood, Metric &metric, Clusterer &clusterer, 
                              const int Nobjects, const double initialEnlargementFactor, const double alpha);
        ~MultiEllipsoidSampler();
        
        virtual void drawWithConstraint(const RefArrayXXd totalSampleOfParameters, const int Nclusters, const RefArrayXi clusterIndices,
                                        const double logWidthInPriorMass, RefArrayXXd drawnSampleOfParameters);
        void computeEllipsoids(const RefArrayXXd totalSampleOfParameters, const int Nclusters, 
                               const RefArrayXi clusterIndices, const double logWidthInPrioMass);
        ArrayXi getNonOverlappingEllipsoidsIndices();
        ArrayXi getOverlappingEllipsoidsIndices();
        vector<Ellipsoid> getEllipsoidsVector();
   

    protected:
      
        void drawFromHyperSphere(Ellipsoid &ellipsoid, RefArrayXd drawnParameters);
        bool intersection(Ellipsoid &ellipsoid1, Ellipsoid &ellipsoid2);
        void findOverlappingEllipsoids();
        bool pointIsInOverlap(Ellipsoid &ellipsoid, const RefArrayXd pointCoordinates);


    private:

        vector<Ellipsoid> ellipsoidsVector;
        ArrayXi NobjectsPerCluster;
        ArrayXi nonOverlappingEllipsoidsIndices;
        ArrayXi overlappingEllipsoidsIndices;
        uniform_real_distribution<> uniform;
        normal_distribution<> normal;           
        int Nobjects;                           // Total number of objects
        int Nellipsoids;                        // Total number of ellipsoids computed
        double initialEnlargementFactor;        // Initial factor for enlargement of ellipsoids
        double alpha;                           // Prior volume shrinkage rate (between 0 and 1)

};

#endif
