// Derived class for sampling from an ellipsoid build
// from a sample of points.
// Created by Enrico Corsaro @ IvS - 28 March 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "HyperQuadricSampler.h"
// Implementation contained in "HyperQuadricSampler.cpp"

#ifndef HYPERQUADRICSAMPLER_H
#define HYPERQUADRICSAMPLER_H

#include <iostream>
#include <iomanip>
#include <random>
#include <functional>
#include <ctime>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "HyperQuadricIntersector.h"
#include "Functions.h"
#include "Metric.h"
#include "Likelihood.h"
#include "Prior.h"


using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXi> RefArrayXi;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class HyperQuadricSampler
{

    public:
       
        HyperQuadricSampler(Prior &prior, Metric &metric, const int Nobjects);
        ~HyperQuadricSampler();
        
        virtual void drawWithConstraint(const int Nclusters, const RefArrayXXd totalSampleOfParameters, RefArrayXi clusterIndices,  // Info from clustering
                                        const double logWidthInPriorMass,                                           // Info from nested sampling
                                        RefArrayXd nestedSampleOfParameters, Likelihood &likelihood) = 0;


    protected:
       
        int Nobjects;
        Prior &prior;
        Metric &metric;


    private:

};

#endif
