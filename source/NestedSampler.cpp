#include "NestedSampler.h"


// NestedSampler::NestedSampler()
//
// PURPOSE: 
//      Constructor. Sets initial information, logEvidence and type 
//      of prior and likelihood distributions to be used. 
//      Increases the number of active nested processes.
//
// INPUT:
//      Nobjects: an integer containing the number of objects to be
//      used in the nested sampling process.
//      ptrPriorsVector: a vector of Prior class objects containing the priors used in the inference.
//      likelihood: a Likelihood class object used for likelihood sampling.
//      metric: a Metric class object to contain the metric used in the problem.
//      clusterer: a Clusterer class object specifying the type of clustering algorithm to be used.
//
// REMARK:
//      The desired model for predictions is to be given initially to 
//      the likelihood object and is not feeded directly inside the 
//      nested sampling process.
//

NestedSampler::NestedSampler(const int Nobjects, vector<Prior*> ptrPriorsVector, Likelihood &likelihood, Metric &metric, Clusterer &clusterer)
: engine(time(0)),
  informationGain(0.0), 
  logEvidence(-DBL_MAX),
  Niterations(0),
  Nobjects(Nobjects),
  ptrPriorsVector(ptrPriorsVector),
  likelihood(likelihood),
  metric(metric),
  clusterer(clusterer)
{
    int totalNdimensions = 0;
    
    for (int i = 0; i < ptrPriorsVector.size(); i++)
    {
        totalNdimensions = ptrPriorsVector[i]->getNdimensions();        // Get the number of dimensions from each type of prior
    }

    Ndimensions = totalNdimensions;
} // END NestedSampler::NestedSampler()








// NestedSampler::~NestedSampler()
//
// PURPOSE: 
//      Destructor. Deletes object of actual Nested sampling computation and
//      decreases the number of active nested processes.
//

NestedSampler::~NestedSampler()
{
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
//      terminationFactor: a double specifying the amount of exceeding information to
//      terminate nested iterations. Default is 1.
//      NiterationsBeforeClustering: number of nested iterations required to recompute
//      the clustering of the points in the sample.
//
// OUTPUT:
//      void
//
// REMARKS: 
//      Eigen Matrices are defaulted column-major. Hence the 
//      nestedSampleOfParameters and posteriorSample are resized as 
//      (Ndimensions, ...), rather than (... , Ndimensions).
//

void NestedSampler::run(const double terminationFactor, const int NiterationsBeforeClustering)
{
    double logWidthInPriorMass;
    double logEvidenceNew;
    double logWeight = 0.0;
    double logMeanEvidence = -DBL_MAX;
    double logMeanLiveEvidence;
    double logTotalLikelihoodOfLivePoints;
    int copy = 0;
    int worst;

    double delta = 1.0/(Nobjects + 1);
    double alpha = delta*Nobjects;
    double beta = 1.0/((Nobjects + 2)*delta);

    // Set up the random number generator. It generates integers random numbers
    // between 0 and Nobjects-1, inclusive. The engine's seed is based on the
    // current time, and a Marsenne Twister pesudo-random generator is used.

    uniform_int_distribution<int> uniform(0, Nobjects-1);


    // Set the sizes of the Eigen Arrays logLikelihood and nestedSampleOfParameters

    logLikelihood.resize(Nobjects);
    nestedSampleOfParameters.resize(Ndimensions, Nobjects);
    ArrayXXd priorSampleOfParameters;
    int NdimensionsPerPrior;               // Number of dimensions having same type of prior
    int actualNdimensions = 0;


    // Initialize all the objects of the process.
    // nestedSampleOfParameters will then contain the sample of coordinates for all the nested objects
    // The first initialization is done by coordinates according to their corresponding prior distribution
    // This process ensures that each object is drawn uniformly from the prior PDF

    for (int i = 0; i < ptrPriorsVector.size(); i++)
    {
        NdimensionsPerPrior = ptrPriorsVector[i]->getNdimensions();
        priorSampleOfParameters.resize(NdimensionsPerPrior, Nobjects);
        ptrPriorsVector[i]->draw(priorSampleOfParameters, Nobjects);
        nestedSampleOfParameters.block(actualNdimensions,0,NdimensionsPerPrior,Nobjects) = priorSampleOfParameters;      
        actualNdimensions += NdimensionsPerPrior;
    }


    // Initialize corresponding likelihood values
    
    ArrayXd nestedSamplePerObject(Ndimensions);
    
    for (int i = 0; i < Nobjects; i++)
    {
        nestedSamplePerObject = nestedSampleOfParameters.col(i);
        logLikelihood(i) = likelihood.logValue(nestedSamplePerObject);
    }


    // Initialize prior mass interval

    logWidthInPriorMass = log(1.0 - exp(-1.0/Nobjects));
    double logTotalWidthInPriorMass = logWidthInPriorMass;
    
    
    // Identify the clusters contained in the initial sample of nested objects

    ArrayXi clusterIndices(Nobjects);
    int Nclusters = clusterer.cluster(nestedSampleOfParameters, clusterIndices);


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
        
        actualLogLikelihoodConstraint = logLikelihood.minCoeff(&worst);
        logWeight = logWidthInPriorMass + actualLogLikelihoodConstraint;                


        // Update the evidence Z and the information Gain
        
        logEvidenceNew = Functions::logExpSum(logEvidence, logWeight);
        informationGain = exp(logWeight - logEvidenceNew) * actualLogLikelihoodConstraint       // Not necessary with new statistical uncertainty
                       + exp(logEvidence - logEvidenceNew) * (informationGain + logEvidence) 
                       - logEvidenceNew;
        logEvidence = logEvidenceNew;

        
        // Save the actual posterior sample and its corresponding likelihood and weights
        
        posteriorSample.col(Niterations) = nestedSampleOfParameters.col(worst);      // save coordinates of worst nested object
        logLikelihoodOfPosteriorSample(Niterations) = actualLogLikelihoodConstraint;       // save corresponding likelihood
        logWeightOfPosteriorSample(Niterations) = logWeight;                         // save corresponding weight ->
                                                                                     // proportional to posterior probability density

        // Compute average likelihood of the actual set of live points
        
        logTotalLikelihoodOfLivePoints = logLikelihood(0);

        for (int m = 1; m < Nobject; m++)
        {
            logTotalLikelihoodOfLivePoints = Functions::logExpSum(totalLogLikelihoodOfLivePoints,logLikelihood(m));
            m++;
        }

        logMeanLikelihoodOfLivePoints = logTotalLikelihoodOfLivePoints - log(Nobjects);
        
        
        // Compute mean evidence and mean live evidence

        logMeanLiveEvidence = logMeanLikelihoodOfLivePoints + Niterations*log(alpha);
        logMeanEvidenceNew = Functions::logExpSum(actualLogLikelihoodConstraint + Niterations*log(alpha) - log(Nobjects), logMeanEvidence);
        logMeanEvidence = logMeanEvidenceNew;
        logMeanTotalEvidence = Functions::logExpSum(logMeanEvidence, logMeanLiveEvidence);
        actualLogTerminationFactor = logMeanLiveEvidence - logMeanEvidence;
        cout << "Skilling's log(Evidence): " << logEvidence << endl;
        cout << "Keeton's log(<Evidence>): " << logMeanEvidence << endl;
        cout << "Keeton's log(<Evidence_Live>): " << logMeanLiveEvidence << endl;
        cout << "Keeton's log(<Total Evidence>): " << logMeanTotalEvidence << endl;
        cout << "------------------------------------" << endl;


        // Replace worst object in favour of a copy of different survivor
        // No replacement if Nobjects == 1. Fundamental step to preserve
        // the randomness of the process

        if (Nobjects > 1)
        {
            do 
            {
                copy = uniform(engine);             // 0 <= copy < Nobjects
            } 
            while (copy == worst);
        }

        nestedSampleOfParameters.col(worst) = nestedSampleOfParameters.col(copy);
        logLikelihood(worst) = logLikelihood(copy);
        

        // If condition is verified, identify the clusters contained in the actual sample of nested objects

        if ((NiterationsBeforeClustering % Niterations) == 0)
            Nclusters = clusterer.cluster(nestedSampleOfParameters, clusterIndices);


        // Evolve the replaced object with the new constraint logLikelihood > actualLogLikelihoodConstraint
        // Compute approximate sampling to find new point verifying the likelihood constraint

        ArrayXXd drawnSampleOfParameters(Ndimensions, 1);
        drawWithConstraint(nestedSampleOfParameters, Nclusters, clusterIndices, logTotalWidthInPriorMass, drawnSampleOfParameters); 
        nestedSampleOfParameters.col(worst) = drawnSampleOfParameters;     // Update new set of parameters in the overall nested sample coordinates
        
        
        // Shrink prior mass interval
        
        logWidthInPriorMass -= 1.0 / Nobjects;

        
        // Update total width in prior mass from beginning to actual nested iteration

        logTotalWidthInPriorMass = Functions::logExpSum(logTotalWidthInPriorMass,logWidthInPriorMass);

        
        // Increase nested loop counter
        
        Niterations++;
        cout << "Niterations: " << Niterations << endl;
    }
    while (Niterations <= 200);
    // while (Niterations <= (terminationFactor * informationGain * Nobjects));   // Termination condition suggested by Skilling 2004
                                                                                  // Run till Niterations >> Nobjects * informationGain
    
 
    // Compute total uncertainty for the Evidence Z
    
    logEvidenceError = sqrt(fabs(informationGain)/Nobjects);

} // END NestedSampler::run()
