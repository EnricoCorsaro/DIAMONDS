#include "NestedSampler.h"

// NestedSampler::NestedSampler()
// PURPOSE: 
//      Class constructor
// INPUT:
//      Ndata = Number of objects for nested sampling
//      Niter = Number of nested iterations
// OUTPUT:
//

NestedSampler::NestedSampler(int Ndata = 100, int Niter = 1000)
: Ndata(Ndata), Niter(Niter)            // What is this for ???
{
    // Set vector sizes
    
    priorM.resize(Ndata);
    param.resize(Ndata);
    logL.resize(Ndata);
    logW.resize(Niter);
    postlogL.resize(Niter);
    postP.resize(Niter);
    results.resize(3);
} // END NestedSampler::NestedSampler




// NestedSampler::run()
// PURPOSE:
//      Start nested sampling computation. Save results in 
//      public data-member vectors "postlogL", "postP", "results"
// INPUT:
// OUTPUT:
//

void NestedSampler::run()
{
    double logwidth;
    double logLstar;
    double H = 0.0;
    double logZ = -DBL_MAX;
    double logZnew;
    int copy;
    int worst;
    MathExtra mathObj;                  // Object providing extra math functions
    
    drawFromPrior();
    logwidth = log(1.0 - exp(-1.0/Ndata));      // initialize prior mass interval

    // Nested sampling loop
    for (int nest = 0; nest < Niter; nest++)
    {
        // Find worst object in the collection
        
        worst = 0;
        for (int i = 1; i < Ndata; i++)
        {
            if (logL.at(i) < logL.at(worst)) worst = i;
        }
        logW.at(worst) = logwidth + logL.at(worst);                
        
        // Update evidence Z and information H
        
        logZnew = mathObj.logExpSum(logZ, logW.at(worst));
        H = updateInformationGain(H, logZ, logZnew, worst);
        logZ = logZnew;

        // Save nested samples for posterior

        postlogL.at(nest) = logL.at(worst);     // save likelihood
        postP.at(nest) = param.at(worst);       // save parameter value
    
        // Replace worst object in favour of a copy of different survivor
        
        srand(time(0));
        do 
        {
            copy = rand() % Ndata;              // 0 <= copy < Ndata
        } 
        while (copy == worst && Ndata > 1);     // do not replace if Ndata = 1

        logLstar = logL.at(worst);
        priorM.at(worst) = priorM.at(copy);
        param.at(worst) = param.at(copy);
        logL.at(worst) = logL.at(copy);
        
        // Evolve the replaced object with the new constraint logL > logLstar
        
        drawFromConstrainedPrior(logLstar, worst);

        // Shrink interval
        
        logwidth -= 1.0 / Ndata;

        // Save the results to public data member
    }
    
    results.at(0) = logZ;
    results.at(1) = sqrt(fabs(H)/Ndata);
    results.at(2) = H;

    return;
} // END NestedSampler::run()




// NestedSampler::drawFromPrior()
// PURPOSE: 
//      Initialize all objects for nested sampling with no constraints on likelihood
// INPUT:
// OUTPUT:
//

void NestedSampler::drawFromPrior()
{
    double param_max, param_min;        // Maximum and minimum parameter values
    double prior_const;                 // 1-dimensional flat prior constant value
    MathExtra mathObj;                  // Object providing extra math functions
            
    param_max = 20;
    param_min = 0;
    prior_const = 1.0 / (param_max - param_min);
    srand(time(0));

    double param0 = 10;                 // Centroid for the Gaussian likelihood
    double sigma = 3.0;                 // Standard deviation for the Gaussian likelihood

    for ( int i = 0; i < Ndata; i++ )
    {
        priorM.at(i) = rand()/(RAND_MAX + 1.);
        param.at(i) = priorM.at(i)/prior_const + param_min;
    }

    mathObj.gaussProfile(param, param0, sigma, 10);	// output stored in mathObj.y
    
    for ( int i = 0; i < Ndata; i++ )
    {
        logL.at(i) = log(mathObj.y.at(i));
    }

    return;
} // END NestedSampler::drawFromPrior()




// NestedSampler::drawFromConstrainedPrior()
// PURPOSE: 
//      Replace worst object in nested sampling with new one subject to 
//      new constraint logL > logL_limit
// INPUT:
//      logL_limit = New constraint on log Likelihood value
//      worst = subscript of edge object in nested iteration
// OUTPUT:
//

void NestedSampler::drawFromConstrainedPrior(double logL_limit, int worst)
{
    double param_max, param_min;        // Maximum and minimum parameter values
    double prior_const;                 // 1-dimensional flat prior constant value
    MathExtra mathObj;                  // Object providing extra math functions
            
    param_max = 20;
    param_min = 0;
    prior_const = 1.0/(param_max - param_min);
    srand(time(0));

    double param0 = 10;                 // Centroid for the Gaussian likelihood
    double sigma = 3.0;                 // Standard deviation for the Gaussian likelihood

    // Find new object subject to constraint logL > LogL_limit

    while (logL.at(worst) < logL_limit)
    {
        priorM.at(worst) = rand()/(RAND_MAX + 1.);
        param.at(worst) = priorM.at(worst)/prior_const + param_min;
        mathObj.gaussProfile(param, param0, sigma, 10);	// output stored in mathObj.y
        logL.at(worst) = log(mathObj.y.at(worst));
    }

    return;

} // END NestedSampler::drawFromConstrainedPrior()




// NestedSampler::updateInformationGain() 
// PURPOSE: 
//      Updates the information gain from old to new evidence
// INPUT:
//      H_old = Old information H
//      logZ_old = old log Evidence
//      logZ_new = new log Evidence
// OUTPUT:
//      New value of information gain H
//

double NestedSampler::updateInformationGain(double H_old, double logZ_old, double logZ_new, int worst)
{    
    return exp(logW[worst] - logZ_new) * logL[worst]
           + exp(logZ_old - logZ_new) * (H_old + logZ_old) - logZ_new;
} // END NestedSampler::updateInformationGain

