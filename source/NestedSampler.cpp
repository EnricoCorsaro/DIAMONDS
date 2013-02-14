#include "NestedSampler.h"



// NestedSampler::NestedSampler()
//
// PURPOSE: 
//      Constructor. Sets initial information, logEvidence and type of prior and likelihood 
//      distributions to be used. Increases the number of active nested processes.
//
// INPUT:
//      variate: a RandomVariate class used as Prior to draw from
//
// OUTPUT:

NestedSampler::NestedSampler(int Ndim)
: , 
  informationH(0.0), 
  logEvidence(-DBL_MAX)
{
    ++nestedCounter;
    cerr << "Nested process initialized" << endl;
    cerr << "Total number of active nested processes: " << nestedCounter << endl;
}






// NestedSampler::~NestedSampler()
//
// PURPOSE: 
//      Destructor. Deletes object of actual Nested sampling computation
//      distributions to be used. Decreases the number of active nested processes.
//
// INPUT:
//
// OUTPUT:

NestedSampler::~NestedSampler()
{
    --nestedCounter;
    cerr << "Nested process deleted" << endl;
    cerr << "Total number of remaining active nested processes: " << nestedCounter << endl;
}







// NestedSampler::getLogEvidence()
//
// PURPOSE:
//      Gets private data member logEvidence
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
//      Gets private data member logEvidenceError
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
//      Gets private data member informationH
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
//      public vectors "logLikelihoodOfPosteriorSample", "posteriorSample"
//
// INPUT:
//      Nobjects : Number of objects for nested sampling
//      Niter : Number of nested iterations
//
// OUTPUT:
//
// REMARKS: Eigen Matrices are defaulted column-major. Hence the parameter and posteriorSample 
//          are resized as (Ndim, ...), rather than (...,Ndim).

void NestedSampler::run(int Nobjects, int Niter)
{
    double logWidthInPriorMass;
    double logLikelihoodConstraint;
    double logEvidenceNew;
    int copy = 0;
    int worst;

    // Set up the random number generator. It generates integers random numbers
    // between 0 and Nobjects-1, inclusive. The engine's seed is based on the
    // current time, and a Marsenne Twister pesudo-random generator is used.
    
    uniform_int_distribution<int> uniform_distribution(0, Nobjects-1);
    mt19937 engine(time(0));
    auto uniform = bind(uniform_distribution, engine);              // Binding uniform distribution to seed value

    // Reset the sizes of the Eigen Arrays
    // Vectors containing objects of a single nested iteration
    parameter.resize(Nobjects);
    logLikelihood.resize(Nobjects);

    // Vectors containing worst objects from all nested iterations
    posteriorSample.resize(Ndim,Niter);
    logLikelihoodOfPosteriorSample.resize(Niter);
    logWeight.resize(Niter);

    // Initialize prior values

    randomVariate.drawNestedValues(parameter, logLikelihood, Nobjects);
    
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

        posteriorSample(nest) = parameter(worst);                          // save parameter value
        logLikelihoodOfPosteriorSample(nest) = logLikelihood(worst);       // save corresponding likelihood

        // Replace worst object in favour of a copy of different survivor
        // No replacement if Nobjects == 1.

        if (Nobjects > 1)
        {
            do 
            {
                copy = uniform();             // 0 <= copy < Nobjects
            } 
            while (copy == worst);
        }

        parameter(worst) = parameter(copy);
        logLikelihood(worst) = logLikelihood(copy);
        
        // Evolve the replaced object with the new constraint logLikelihood > logLikelihoodConstraint
        
        randomVariate.drawNestedValueWithConstraint(parameter(worst), logLikelihood(worst), logLikelihoodConstraint);

        // Shrink interval

        logWidthInPriorMass -= 1.0 / Nobjects;
    }
    
    // Compute uncertainty on the log of the Evidence Z

    logEvidenceError = sqrt(fabs(informationH)/Nobjects);

    return;
}
