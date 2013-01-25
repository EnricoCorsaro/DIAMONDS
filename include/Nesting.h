// Nesting class for setting up (mono)nested sampling analysis.
// This file presents the public interface of the class.
// Member-function definitions appear in Nesting.cpp.
// Created by Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "Nesting.h"

#include <cfloat>
#include "MathExtra.h"
#include "FileProcess.h"


class Nesting
{

    public:
        
        // Data members accessible from outside
        vector<double> param;           // parameter values (the free parameters of the problem)
        vector<double> postlogL;        // likelihood samples from nested sampling
        vector<double> postP;           // parameter values corresponding to posterior sample
        vector<double> results;         // output logZ, logZ_err, information H

        // Constructor
        Nesting ( int nobj, int nnest );


	private:

        // Private data members
        int n_obj;
        int n_nest;
        vector<double> priorM;          // prior mass values between 0 and last nested boundary
        vector<double> logL;            // log-likelihood values
        vector<double> logW;            // weight = width * likelihood

        // Core of the nested sampling algorithm
        void nestedSampling ();

        // Initialize all objects for nested sampling with no constraints on likelihood
        void priorSampling ();
        
        // Replace worst object in nested sampling with new one subject to constraint L > L*
        void constrainSampling ( double logL_limit, int worst );
        
        // Computes the information gain from old to new evidence
        double informationGain ( double logZ_old, double logZ_new, int worst );
}   ;
