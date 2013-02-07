#include "NestedSampler.h"

// NestedSampler::NestedSampler()
// PURPOSE: 
//      Constructor. Sets initial information, logEvidence and type of prior and likelihood 
//      distributions to be used.
// INPUT:
//      variate: a RandomVariate class used as Prior to draw from
// OUTPUT:

NestedSampler::NestedSampler(RandomVariate &variate)
: randomVariate(variate), informationH(0.0), logEvidence(-DBL_MAX)
{

}




// NestedSampler::getLogEvidence()
// PURPOSE:
//      Gets private data member logEvidence
// INPUT:
// OUTPUT:

double NestedSampler::getLogEvidence()
{
    return logEvidence;
}




// NestedSampler::getLogEvidenceError()
// PURPOSE:
//      Gets private data member logEvidenceError
// INPUT:
// OUTPUT:

double NestedSampler::getLogEvidenceError()
{
    return logEvidenceError;
}




// NestedSampler::getInformationH()
// PURPOSE:
//      Gets private data member informationH
// INPUT:
// OUTPUT:

double NestedSampler::getInformationH()
{
    return informationH;
}




// NestedSampler::run()
// PURPOSE:
//      Start nested sampling computation. Save results in 
//      public vectors "logLikelihoodOfPosteriorSample", "posteriorSample"
// INPUT:
//      Nobjects = Number of objects for nested sampling
//      Niter = Number of nested iterations
// OUTPUT:

void NestedSampler::run(int Nobjects, int Niter)
{
    double logWidthInPriorMass;
    double logLikelihoodConstraint;
    double logEvidenceNew;
    int copy;
    int worst;

    // Set up the random number generator. It generates random numbers
    // between 0 and Nobjects-1, inclusive. The engine's seed is based on the 
    // current time and uses a Marsenne Twister pesudo-random generator.
    uniform_int_distribution<int> uniform_distribution(0, Nobjects-1);
    mt19937 engine(time(0));
    auto uniform = bind(uniform_distribution, engine);              // Binding uniform distribution to seed value

    // Set vector sizes
    parameter.resize(Nobjects);
    logLikelihood.resize(Nobjects);
    logWeight.resize(Niter);
    posteriorSample.resize(Niter);
    logLikelihoodOfPosteriorSample.resize(Niter);

    // Initialize prior values
    randomVariate.drawNestedValues(parameter, logLikelihood, Nobjects);
    
    // Initialize prior mass interval
    logWidthInPriorMass = log(1.0 - exp(-1.0/Nobjects));  

    // Nested sampling loop
    for (int nest = 0; nest < Niter; nest++)
    {
        // Find worst object in the collection
        worst = 0;
        for (int i = 1; i < Nobjects; i++)
        {
            if (logLikelihood.at(i) < logLikelihood.at(worst))
            {
                worst = i;
            }
        }
        logWeight.at(worst) = logWidthInPriorMass + logLikelihood.at(worst);                
        
        // Update evidence Z and information H
        logEvidenceNew = MathExtra::logExpSum(logEvidence, logWeight.at(worst));
        informationH = updateInformationGain(informationH, logEvidence, logEvidenceNew, worst);
        logEvidence = logEvidenceNew;

        // Save the posterior sample and its corresponding likelihood
        posteriorSample.at(nest) = parameter.at(worst);                         // save parameter value
        logLikelihoodOfPosteriorSample.at(nest) = logLikelihood.at(worst);  // save corresponding likelihood
    
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
        cout << "Here" << endl;

        logLikelihoodConstraint = logLikelihood.at(worst);
        parameter.at(worst) = parameter.at(copy);
        logLikelihood.at(worst) = logLikelihood.at(copy);
        
        // Evolve the replaced object with the new constraint logLikelihood > logLikelihoodConstraint
        randomVariate.drawNestedValueWithConstraint(parameter.at(worst), logLikelihood.at(worst), logLikelihoodConstraint);
        
        // Shrink interval
        logWidthInPriorMass -= 1.0 / Nobjects;
    }
    
    // Compute uncertainty on the log of the Evidence Z
    logEvidenceError = sqrt(fabs(informationH)/Nobjects);

    return;
}




// NestedSampler::updateInformationGain() 
// PURPOSE: 
//      Updates the information gain from old to new evidence
// INPUT:
//      H_old = Old information H
//      logEvidence_old = old log Evidence
//      logEvidence_new = new log Evidence
// OUTPUT:
//      New value of information gain H

double NestedSampler::updateInformationGain(double H_old, double logEvidence_old, double logEvidence_new, int worst)
{    
    return exp(logWeight[worst] - logEvidence_new) * logLikelihood[worst]
           + exp(logEvidence_old - logEvidence_new) * (H_old + logEvidence_old) - logEvidence_new;
}
