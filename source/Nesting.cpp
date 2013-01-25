// Nesting member-function definitions. This file contains
// implementations of the member functions prototyped in GradeBook.h.
// Created by Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Source file "Nesting.cpp"
#include "MathExtra.h"
#include "Nesting.h"

        // Constructor
        Nesting::Nesting ( int nobj = 100, int nnest = 1000 )
        {
            n_obj   = nobj;             // Number of objects for nested sampling
            n_nest  = nnest;            // Number of nested iterations
            
            // Set vector sizes
            priorM.resize(n_obj);
            param.resize(n_obj);
            logL.resize(n_obj);
            logW.resize(n_nest);
            postlogL.resize(n_nest);
            postP.resize(n_nest);
            results.resize(3);

            // Start computation
            nestedSampling();
        }

        // Core of the nested sampling algorithm
        void Nesting::nestedSampling ()
        {
            double logwidth;
            double logLstar;
            double H = 0.0;
            double logZ = -DBL_MAX;
            double logZnew;
            int copy;
            int worst;
	        MathExtra mathObj;                  // Object providing extra math functions
            
            priorSampling();
            logwidth = log(1.0 - exp(-1.0/n_obj));      // initialize prior mass interval

            // Nested sampling loop
            for ( int nest = 0; nest < n_nest; nest++ )
            {

                // Find worst object in the collection
                worst = 0;
                for ( int i = 1; i < n_obj; i++ )
                    {
                        if ( logL.at(i) < logL.at(worst) )
                            worst = i;
                    }
                logW.at(worst) = logwidth + logL.at(worst);                
                
                // Update evidence Z and information H
                logZnew = mathObj.logExpSum ( logZ, logW.at(worst) );
                H = informationGain( logZ, logZnew, worst );
                logZ = logZnew;

                // Save nested samples for posterior
                postlogL.at(nest) = logL.at(worst);        // save likelihood
                postP.at(nest) = param.at(worst);       // save parameter value
            
                // Replace worst object in favour of a copy of different survivor
                srand(time(0));
                do 
                {
                    copy = rand() % n_obj;              // 0 <= copy < n_obj
                } while ( copy == worst && n_obj > 1 ); // do not replace if n_obj = 1

                logLstar = logL.at(worst);
                priorM.at(worst) = priorM.at(copy);
                param.at(worst) = param.at(copy);
                logL.at(worst) = logL.at(copy);
                
                // Evolve the replaced object with the new constraint logL > logLstar
                constrainSampling ( logLstar, worst );

                // Shrink interval
                logwidth -= 1.0 / n_obj;

                // Save the results to public data member
            }
            
            results.at(0) = logZ;
            results.at(1) = sqrt(fabs(H)/n_obj);
            results.at(2) = H;

            return;
        }

        // Initialize all objects for nested sampling with no constraints on likelihood
        void Nesting::priorSampling ()
        {
	        double param_max, param_min;        // Maximum and minimum parameter values
	        double prior_const;                 // 1-dimensional flat prior constant value
	        MathExtra mathObj;                  // Object providing extra math functions
            	    
            // BEGIN - Example of flat prior and a Gaussian likelihood
	        param_max = 20;
	        param_min = 0;
	        prior_const = pow( (param_max - param_min), -1. );
	        srand(time(0));
	
	        double param0 = 10;                 // Centroid for the Gaussian likelihood
	        double sigma = 3.0;                 // Standard deviation for the Gaussian likelihood

	        for ( int i = 0; i < n_obj; i++ )
	        {
		        priorM.at(i) = rand()/(RAND_MAX + 1.);
		        param.at(i) = priorM.at(i)/prior_const + param_min;
		        //cout << x.at(i) << endl;
	        }

	        mathObj.gaussProfile ( param, param0, sigma, 10 );	// output stored in mathObj.y
	        
            for ( int i = 0; i < n_obj; i++ )
                logL.at(i) = log(mathObj.y.at(i));
	        // END - Example of flat prior and a Gaussian likelihood

            return;
        }
        
        // Replace worst object in nested sampling with new one subject to constraint L > L*
        void Nesting::constrainSampling ( double logL_limit, int worst )
        {
	        double param_max, param_min;        // Maximum and minimum parameter values
	        double prior_const;                 // 1-dimensional flat prior constant value
	        MathExtra mathObj;                  // Object providing extra math functions
            	    
            // BEGIN - Example of flat prior and a Gaussian likelihood
	        param_max = 20;
	        param_min = 0;
	        prior_const = pow( (param_max - param_min), -1. );
	        srand(time(0));
	
	        double param0 = 10;                 // Centroid for the Gaussian likelihood
	        double sigma = 3.0;                 // Standard deviation for the Gaussian likelihood

            // Find new object subject to constraint logL > LogL_limit
            while (logL.at(worst) < logL_limit)
            {
		        priorM.at(worst) = rand()/(RAND_MAX + 1.);
		        param.at(worst) = priorM.at(worst)/prior_const + param_min;
		        //cout << x.at(i) << endl;

	            mathObj.gaussProfile ( param, param0, sigma, 10 );	// output stored in mathObj.y
                logL.at(worst) = log(mathObj.y.at(worst));
            }
	        // END - Example of flat prior and a Gaussian likelihood

            return;
        }
        
        // Computes the information gain from old to new evidence
        double Nesting::informationGain ( double logZ_old, double logZ_new, int worst )
        {
            double info;
            info = exp ( logW.at(worst) - logZ_new ) * logL.at(worst)
                + exp ( logZ_old - logZ_new ) * ( info + logZ_old ) - logZ_new;

            return info;
        }
}   ;
