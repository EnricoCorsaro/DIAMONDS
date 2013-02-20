#include "NestedSampler.h"


// NestedSampler::NestedSampler()
//
// PURPOSE: 
//      Constructor. Sets initial information, logEvidence and type 
//      of prior and likelihood distributions to be used. 
//      Increases the number of active nested processes.
//
// INPUT:
//      prior: a prior class object used as Prior to draw from
//      likelihood: a likelihood class object used for likelihood sampling
//
// REMARK:
//      The desired model for predictions is to be given initially to 
//      the likelihood object and is not used directly inside the 
//      nested sampling process.
//

NestedSampler::NestedSampler(Prior &prior, Likelihood &likelihood)
: informationGain(0.0), 
  logEvidence(-DBL_MAX),
  nestIteration(0),
  prior(prior),
  likelihood(likelihood)
{
    // ??? Not working
//    NestedSampler::nestedCounter++;
//    cerr << "Nested process initialized" << endl;
//    cerr << "Total number of active nested processes: " << NestedSampler::nestedCounter << endl;
   
} // END NestedSampler::NestedSampler()






// NestedSampler::~NestedSampler()
//
// PURPOSE: 
//      Destructor. Deletes object of actual Nested sampling computation and
//      decreases the number of active nested processes.
//

NestedSampler::~NestedSampler()
{
    // ??? Not working
//    NestedSampler::nestedCounter--;
//    cerr << "Nested process deleted" << endl;
//    cerr << "Total number of remaining active nested processes: " << NestedSampler::nestedCounter << endl;
} // END NestedSampler::~NestedSampler()







// NestedSampler::getLogEvidence()
//
// PURPOSE:
//      Get private data member logEvidence.
//
// OUTPUT:
//      A double containing the natural logarithm of the final evidence.
//

double NestedSampler::getLogEvidence()
{
    return logEvidence;
} // END NestedSampler::logEvidence()







// NestedSampler::getLogEvidenceError()
//
// PURPOSE:
//      Get private data member logEvidenceError.
//
// OUTPUT:
//      A double containing the natural logarithm of the 
//      final error on the evidence.
//

double NestedSampler::getLogEvidenceError()
{
    return logEvidenceError;
} // END NestedSampler::logEvidenceError()








// NestedSampler::getInformationGain()
//
// PURPOSE:
//      Get private data member informationGain.
//
// OUTPUT:
//      A double containing the final amount of
//      information gain in moving from prior to posterior.
//
double NestedSampler::getInformationGain()
{
    return informationGain;
} // END NestedSampler::getInformationGain()










// NestedSampler::getNestIteration()
//
// PURPOSE:
//      Get private data member nestIteration.
//
// OUTPUT:
//      An integer containing the final number of
//      nested loop iterations.
//
int NestedSampler::getNestIteration()
{
    return nestIteration;
} // END NestedSampler::getNestIteration()











// NestedSampler::run()
//
// PURPOSE:
//      Start nested sampling computation. Save results in 
//      arrays "logLikelihoodOfPosteriorSample", "posteriorSample".
//
// OUTPUT:
//      void
//
// REMARKS: 
//      Eigen Matrices are defaulted column-major. Hence the 
//      nestedParameters and posteriorSample are resized as 
//      (Ndim , ...), rather than (... , Ndim).
//

void NestedSampler::run()
{
    double logWidthInPriorMass;
    double logLikelihoodConstraint;
    double logEvidenceNew;
    double logWeight = 0;
    double exceedFactor = 1.2;                  // Defines the termination condition for the nested sampling loop !!! Very critical !!!
    int Nobjects = prior.getNobjects();
    int Ndimensions = prior.getNdimensions();
    int copy = 0;
    int worst;

    // Set up the random number generator. It generates integers random numbers
    // between 0 and Nobjects-1, inclusive. The engine's seed is based on the
    // current time, and a Marsenne Twister pesudo-random generator is used.
    uniform_int_distribution<int> uniform_distribution(0, Nobjects-1);
    mt19937 engine(time(0));
    auto uniform = bind(uniform_distribution, engine);              // Binding uniform distribution to seed value

    // Set the sizes of the Eigen Arrays logLikelihood and nestedParameters
    logLikelihood.resize(Nobjects);
    nestedParameters.resize(Ndimensions, Nobjects);

    // Initialize the objects
    prior.draw(nestedParameters);           // nestedParameters will then contain the sample of parameters for nested sampling
    
    // Initialize corresponding likelihood values
    ArrayXd objectParameters(Ndimensions);
    
    for (int i = 0; i < Nobjects; i++)
    {
        objectParameters = nestedParameters.col(i);
        logLikelihood(i) = likelihood.logValue(objectParameters);
    }

    // Initialize prior mass interval
    logWidthInPriorMass = log(1.0 - exp(-1.0/Nobjects));

    // Nested sampling loop
    do 
    {
        // Resizing array dimensions to the actual number of nested iterations
        // conservativeResize allows dinamic resizing of Eigen Arrays, while keeping the previous values untouched
        posteriorSample.conservativeResize(Ndimensions, nestIteration + 1);        // conservative resize to column number only
        logLikelihoodOfPosteriorSample.conservativeResize(nestIteration + 1);
        
        // Find worst object in the collection. The likelihood of this object
        // defines the constraint when drawing new objects later on.
        logLikelihoodConstraint = logLikelihood.minCoeff(&worst);
        logWeight = logWidthInPriorMass + logLikelihoodConstraint;                
        
        // Update the evidence Z and the information Gain
        logEvidenceNew = MathExtra::logExpSum(logEvidence, logWeight);
        informationGain = exp(logWeight - logEvidenceNew) * logLikelihoodConstraint
                       + exp(logEvidence - logEvidenceNew) * (informationGain + logEvidence) 
                       - logEvidenceNew;
        logEvidence = logEvidenceNew;

        // Set actual array size and save the posterior sample and its corresponding likelihood
        posteriorSample.col(nestIteration) = nestedParameters.col(worst);              // save parameter value
        logLikelihoodOfPosteriorSample(nestIteration) = logLikelihoodConstraint;       // save corresponding likelihood

        // Replace worst object in favour of a copy of different survivor
        // No replacement if Nobjects == 1. Fundamental step to preserve
        // the randomness of the process
        if (Nobjects > 1)
        {
            do 
            {
                copy = uniform();             // 0 <= copy < Nobjects
            } 
            while (copy == worst);
        }

        nestedParameters.col(worst) = nestedParameters.col(copy);
        logLikelihood(worst) = logLikelihood(copy);
        
        // Evolve the replaced object with the new constraint logLikelihood > logLikelihoodConstraint
        objectParameters = nestedParameters.col(worst);
        prior.drawWithConstraint(objectParameters, likelihood);     // Array objectParameters is updated after call to function
        nestedParameters.col(worst) = objectParameters;             // Update new set of parameters in nestedParameters array
        
        // Shrink interval
        logWidthInPriorMass -= 1.0 / Nobjects;

        // Increase nested loop counter
        nestIteration++;
        cout << "nestIteration: " << nestIteration << endl;
        cout << "Information Gain * Nobjects : " << informationGain * Nobjects << endl;
    }
    while (nestIteration <= (exceedFactor * informationGain * Nobjects));   // Termination condition suggested by Skilling 2004
                                                                            // Run till nestIteration >> Nobjects * informationGain
    // Compute uncertainty on the log of the Evidence Z
    logEvidenceError = sqrt(fabs(informationGain)/Nobjects);

} // END NestedSampler::run()
