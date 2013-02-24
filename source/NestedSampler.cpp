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
  Niterations(0),
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










// NestedSampler::getNiterations()
//
// PURPOSE:
//      Get private data member Niterations.
//
// OUTPUT:
//      An integer containing the final number of
//      nested loop iterations.
//
int NestedSampler::getNiterations()
{
    return Niterations;
} // END NestedSampler::getNiterations()










// NestedSampler::run()
//
// PURPOSE:
//      Start nested sampling computation. Save results in Eigen
//      Arrays logLikelihoodOfPosteriorSample, posteriorSample,
//      logWeightOfPosteriorSample.
//
// INPUT:
//      Nobjects: an integer containing the number of objects to be
//      used in the nested sampling process.
//
// OUTPUT:
//      void
//
// REMARKS: 
//      Eigen Matrices are defaulted column-major. Hence the 
//      nestedSampleOfParameters and posteriorSample are resized as 
//      (Ndim , ...), rather than (... , Ndim).
//

void NestedSampler::run(const int Nobjects)
{
    double logWidthInPriorMass;
    double logLikelihoodConstraint;
    double logEvidenceNew;
    double logWeight = 0.0;
    double exceedFactor = 1.2;                  // Defines the termination condition for the nested sampling loop 
                                                // !!! Very critical !!!
    int Ndimensions = prior.getNdimensions();
    int copy = 0;
    int worst;


    // Set up the random number generator. It generates integers random numbers
    // between 0 and Nobjects-1, inclusive. The engine's seed is based on the
    // current time, and a Marsenne Twister pesudo-random generator is used.

    uniform_int_distribution<int> uniform_distribution(0, Nobjects-1);
    mt19937 engine(time(0));
    auto uniform = bind(uniform_distribution, engine);              // Binding uniform distribution to seed value


    // Set the sizes of the Eigen Arrays logLikelihood and nestedSampleOfParameters

    logLikelihood.resize(Nobjects);
    nestedSampleOfParameters.resize(Ndimensions, Nobjects);


    // Initialize the objects
    // nestedSampleOfParameters will then contain the sample of parameters for nested sampling

    prior.draw(nestedSampleOfParameters, Nobjects);         


    // Initialize corresponding likelihood values
    
    ArrayXd nestedSamplePerObject(Ndimensions);
    
    for (int i = 0; i < Nobjects; i++)
    {
        nestedSamplePerObject = nestedSampleOfParameters.col(i);
        logLikelihood(i) = likelihood.logValue(nestedSamplePerObject);
    }


    // Initialize prior mass interval

    logWidthInPriorMass = log(1.0 - exp(-1.0/Nobjects));


    // Nested sampling loop

    do 
    {
        // Resizing array dimensions to the actual number of nested iterations
        // conservativeResize allows dinamic resizing of Eigen Arrays, while keeping the previous values untouched
    
        posteriorSample.conservativeResize(Ndimensions, Niterations + 1);        // conservative resize to column number only
        logLikelihoodOfPosteriorSample.conservativeResize(Niterations + 1);
        logWeightOfPosteriorSample.conservativeResize(Niterations + 1);
        

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

        
        // Save the actual posterior sample and its corresponding likelihood and weights
        
        posteriorSample.col(Niterations) = nestedSampleOfParameters.col(worst);              // save parameter value
        logLikelihoodOfPosteriorSample(Niterations) = logLikelihoodConstraint;       // save corresponding likelihood
        logWeightOfPosteriorSample(Niterations) = logWeight;                         // save corresponding weight ->
                                                                                     // proportional to posterior probability density

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

        nestedSampleOfParameters.col(worst) = nestedSampleOfParameters.col(copy);
        logLikelihood(worst) = logLikelihood(copy);
        

        // Evolve the replaced object with the new constraint logLikelihood > logLikelihoodConstraint
        
        nestedSamplePerObject = nestedSampleOfParameters.col(worst);
        prior.drawWithConstraint(nestedSamplePerObject, likelihood);     // Array nestedSamplePerObject is updated after call to function
        nestedSampleOfParameters.col(worst) = nestedSamplePerObject;     // Update new set of parameters in nestedSampleOfParameters array
        

        // Shrink interval
        
        logWidthInPriorMass -= 1.0 / Nobjects;

        
        // Increase nested loop counter
        
        Niterations++;
    }
    while (Niterations <= 300);   // Termination condition suggested by Skilling 2004
    // while (Niterations <= (exceedFactor * informationGain * Nobjects));   // Termination condition suggested by Skilling 2004
                                                                            // Run till Niterations >> Nobjects * informationGain
    
    cout << "Niterations: " << Niterations << endl;
    cout << "Information Gain * Nobjects : " << informationGain * Nobjects << endl;
 
    // Compute uncertainty on the log of the Evidence Z
    
    logEvidenceError = sqrt(fabs(informationGain)/Nobjects);

} // END NestedSampler::run()
