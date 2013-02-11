#include "NestedSampler.h"

// NestedSampler::NestedSampler()
//
// PURPOSE: 
//      Class constructor
//
// INPUT:
//      variate: a RandomVariate class used as Prior to draw from
//
// OUTPUT:

NestedSampler::NestedSampler(RandomVariate &variate)
: randomVariate(variate), informationH(0.0), logEvidence(-DBL_MAX)
{

}






// NestedSampler::getLogEvidence()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

double NestedSampler::getLogEvidence()
{
    return logEvidence;
}








// NestedSampler::getLogEvidenceError()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

double NestedSampler::getLogEvidenceError()
{
    return logEvidenceError;
}







// NestedSampler::getInformationH()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

double NestedSampler::getInformationH()
{
    return informationH;
}








// NestedSampler::run()
//
// PURPOSE:
//      Start nested sampling computation. Save results in 
//      public vectors "logLikelihoodOfPosteriorSample", "posteriorSample", "results"
//
// INPUT:
//      Nobjects = Number of objects for nested sampling
//      Niter = Number of nested iterations
//
// OUTPUT:
//
// REMARKS: Eigen Matrices are defaulted column-major. Hence the param and posteriorSample 
//          are resized as (Ndim, ...), rather than (...,Ndim).
//          FIXME: we're using integers for indexing the arrays. Size of arrays may
//                 be too large for this. Check out std::ptrdiff_t.

void NestedSampler::run(int Nobjects, int Niter)
{
    double logWidthInPriorMass;
    double logLikelihoodConstraint;
    double logEvidenceNew;
    int copy=0;
    int worst;

    // Set up the random number generator. It generates integers random numbers
    // between 0 and Nobjects-1, inclusive. The engine's seed is based on the 
    // current time.
    
    uniform_int_distribution<int> uniform_distribution(0, Nobjects-1);
    mt19937 engine(time(0));
    auto uniform = bind(uniform_distribution, engine);

    // Reset the sizes of the Eigen Arrays
    
    param.resize(randomVariate.getNdim(), Nobjects);
    logLikelihood.resize(Nobjects);
    logWeight.resize(Niter);
    posteriorSample.resize(randomVariate.getNdim(), Niter);
    logLikelihoodOfPosteriorSample.resize(Niter);

    // Initialize prior values
    
    randomVariate.drawNestedValues(param, logLikelihood, Nobjects);
    
    // Initialize prior mass interval
    
    logWidthInPriorMass = log(1.0 - exp(-1.0/Nobjects));  

    // Nested sampling loop
        
    for (int nest = 0; nest < Niter; nest++)
    {
        // Find worst object in the collection. The likelihood of this object
        // defines the constraint when drawing new objects later on.
        
        logLikelihoodConstraint = logLikelihood.minCoeff(&worst);
        logWeight(worst) = logWidthInPriorMass + logLikelihood(worst);                
        
        // Update the evidence Z and the information H
        
        logEvidenceNew = MathExtra::logExpSum(logEvidence, logWeight(worst));
        informationH = exp(logWeight(worst) - logEvidenceNew) * logLikelihood(worst)
                       + exp(logEvidence - logEvidenceNew) * (informationH + logEvidence) 
                       - logEvidenceNew;
        logEvidence = logEvidenceNew;

        // Save the posterior sample and its corresponding likelihood

        posteriorSample.col(nest) = param.col(worst);                       // save parameter value
        logLikelihoodOfPosteriorSample(nest) = logLikelihood(worst);        // save corresponding likelihood
    
        // Replace worst object in favour of a copy of different survivor
        // No replacement if Nobjects == 1.
        
        if (Nobjects > 1)
        {
            do 
            {
                copy = uniform();              // 0 <= copy < Nobjects
            } 
            while (copy == worst);
        }

        param.col(worst) = param.col(copy);
        logLikelihood(worst) = logLikelihood(copy);
        
        // Evolve the replaced object with the new constraint logLikelihood > logLikelihoodConstraint
        
        randomVariate.drawNestedValueWithConstraint(param.col(worst), logLikelihood(worst), logLikelihoodConstraint);
        
        // Shrink interval
        
        logWidthInPriorMass -= 1.0 / Nobjects;

        // Save the results to public data member
    }
    
    // Compute uncertainty on the log of the Evidence Z
    
    logEvidenceError = sqrt(fabs(informationH)/Nobjects);

    return;
}


