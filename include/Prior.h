// Class for setting and drawing from prior distribution
// Enrico Corsaro @ IvS - 13 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "Prior.h"
// Implementation contained in "Prior.cpp"

#ifndef PRIOR_H
#define PRIOR_H

#include <random>
#include <Eigen/Core>
#include "MathExtra.h"
#include "Model.h"
#include "Likelihood.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;

class Prior
{

    public:

        Prior(int Ndim);
        int getNdim();
        
        // Functions for defining parameter intervals
        void setBoundaries(const RefArrayXd min, const RefArrayXd max); 
        double getMinimum();
        double getMaximum();

        // Functions for uniform prior distribution (flat, uninformative)
        void uniformPrior();
        double getUniformFactor();
        void drawFromUniformPrior(RefArrayXXd NestedSampleParameters, const double Nobjects);
        void drawFromUniformPriorWithConstraint(RefArrayXd NestedSampleParameter, const double logLikelihoodConstraint);


    private:

        uniform_real_distribution<> uniform;
        mt19937 engine;
        int Ndim;
        double uniformFactor;
        ArrayXd minimum;
        ArrayXd maximum;

}; // END class Prior

#endif
