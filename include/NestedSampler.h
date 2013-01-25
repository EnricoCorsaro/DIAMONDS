
#ifndef NESTEDSAMPLER_H
#define NESTEDSAMPLER_H

#include <cfloat>
#include "MathExtra.h"
#include "FileProcess.h"


class NestedSampler
{
    public:
        
        NestedSampler(int Ndata, int Niter);    // constructor
        void run();
        vector<double> param;                   // parameter values (the free parameters of the problem)
        vector<double> postlogL;                // likelihood samples from nested sampling
        vector<double> postP;                   // parameter values corresponding to posterior sample
        vector<double> results;                 // output logZ, logZ_err, information H

	private:

        int Ndata;
        int Niter;
        
        vector<double> priorM;          // prior mass values between 0 and last nested boundary
        vector<double> logL;            // log-likelihood values
        vector<double> logW;            // weight = width * likelihood

        void drawFromPrior();
        void drawFromConstrainedPrior(double logL_limit, int worst);
        double informationGain(double oldH, double logZ_old, double logZ_new, int worst);
};

#endif
