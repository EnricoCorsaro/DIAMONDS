#include "NestedSampler.h"


// NestedSampler::NestedSampler()
//
// PURPOSE: 
//      Constructor. Sets initial information, logEvidence and type 
//      of prior and likelihood distributions to be used. 
//      Increases the number of active nested processes.
//
// INPUT:
//      printOnTheScreen: Boolean value specifying whether the results are to 
//                         be printed on the screen or not.
//      Nobjects:         Integer containing the number of objects to be
//                         used in the nested sampling process.
//      ptrPriors:        Vector of pointers to Prior class objects
//      likelihood:       Likelihood class object used for likelihood sampling.
//      metric:           Metric class object to contain the metric used in the problem.
//      clusterer:        Clusterer class object specifying the type of clustering algorithm to be used.
//
// REMARK:
//      The desired model for predictions is to be given initially to 
//      the likelihood object and is not feeded directly inside the 
//      nested sampling process.
//

NestedSampler::NestedSampler(const bool printOnTheScreen, const int Nobjects, vector<Prior*> ptrPriors, 
                             Likelihood &likelihood, Metric &metric, Clusterer &clusterer)
: ptrPriors(ptrPriors),
  likelihood(likelihood),
  metric(metric),
  clusterer(clusterer),
  printOnTheScreen(printOnTheScreen),
  Nobjects(Nobjects),
  engine(time(0)),
  uniform(0.0, 1.0),
  Niterations(0),
  informationGain(0.0), 
  logEvidence(-DBL_MAX),
  logMeanEvidence(-DBL_MAX)

{
    int totalNdimensions = 0;
    
    for (int i = 0; i < ptrPriors.size(); i++)
    {
        // Get the number of dimensions from each type of prior

        totalNdimensions = ptrPriors[i]->getNdimensions(); 
    }

    Ndimensions = totalNdimensions;
    constant1 = 1.0/(Nobjects + 1);
    constant2 = constant1*Nobjects;
    constant3 = 1.0/((Nobjects + 2)*constant1);

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
//      maxNdrawAttempts: the maximum number of attempts allowed when drawing from a single ellipsoid.
//
// OUTPUT:
//      void
//
// REMARKS: 
//      Eigen Matrices are defaulted column-major. Hence the 
//      nestedSampleOfParameters and posteriorSample are resized as 
//      (Ndimensions, ...), rather than (... , Ndimensions).
//

void NestedSampler::run(const double terminationFactor, const int NiterationsBeforeClustering, const int maxNdrawAttempts)
{
    int copy = 0;
    int worst;
    int startTime = time(0);
    int NdimensionsPerPrior;               // Number of dimensions having same type of prior
    int actualNdimensions = 0;
    double logWidthInPriorMass;
    double logEvidenceNew;
    double logWeight = 0.0;
    double logMeanEvidenceNew;
    double logMeanLiveEvidence;
    double logTotalLikelihoodOfLivePoints = -DBL_MAX;
    double actualTerminationFactor;


    // Set up the random number generator. It generates integers random numbers
    // between 0 and Nobjects-1, inclusive. The engine's seed is based on the
    // current time, and a Marsenne Twister pesudo-random generator is used.

    uniform_int_distribution<int> uniform(0, Nobjects-1);


    // Set the sizes of the Eigen Arrays logLikelihood and nestedSampleOfParameters

    logLikelihood.resize(Nobjects);
    nestedSampleOfParameters.resize(Ndimensions, Nobjects);
    ArrayXXd priorSampleOfParameters;


    // Initialize all the objects of the process.
    // nestedSampleOfParameters will then contain the sample of coordinates for all the nested objects
    // The first initialization is done by coordinates according to their corresponding prior distribution
    // This process ensures that each object is drawn uniformly from the prior PDF

    for (int i = 0; i < ptrPriors.size(); i++)
    {
        NdimensionsPerPrior = ptrPriors[i]->getNdimensions();
        priorSampleOfParameters.resize(NdimensionsPerPrior, Nobjects);
        ptrPriors[i]->draw(priorSampleOfParameters);
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

    int Nclusters;
    ArrayXi clusterIndices(Nobjects);


    // Nested sampling loop

    do 
    {
        if (printOnTheScreen && (Niterations !=0))
        {
            for (int l=0; l < 20; l++)
            {
                cerr << endl;
            }
        }


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
        logLikelihoodOfPosteriorSample(Niterations) = actualLogLikelihoodConstraint; // save corresponding likelihood
        logWeightOfPosteriorSample(Niterations) = logWeight;                         // save corresponding weight...
                                                                                     // ...proportional to posterior probability density

        // Compute average likelihood of the actual set of live points
        
        logTotalLikelihoodOfLivePoints = logLikelihood(0);

        for (int m = 1; m < Nobjects; m++)
        {
            logTotalLikelihoodOfLivePoints = Functions::logExpSum(logTotalLikelihoodOfLivePoints,logLikelihood(m));
        }

        logMeanLikelihoodOfLivePoints = logTotalLikelihoodOfLivePoints - log(Nobjects);
        
        
        // Compute mean evidence and mean live evidence

        logMeanLiveEvidence = logMeanLikelihoodOfLivePoints + Niterations*log(constant2);
        logMeanEvidenceNew = Functions::logExpSum(actualLogLikelihoodConstraint + Niterations*log(constant2) - log(Nobjects), logMeanEvidence);
        logMeanEvidence = logMeanEvidenceNew;
        logMeanTotalEvidence = Functions::logExpSum(logMeanEvidence, logMeanLiveEvidence);


        // Evaluate termination condition
        
        actualTerminationFactor = exp(logMeanLiveEvidence - logMeanEvidence);
        
        
        // Print current information on the screen, if required

        if (printOnTheScreen)
        {
            cerr << "=========================================" << endl;
            cerr << "Information on Nesting process" << endl;
            cerr << "=========================================" << endl;
            cerr << "Niterations: " << Niterations << endl;
            cerr << "Total Width In Prior Mass: " << exp(logTotalWidthInPriorMass) << endl;
            cerr << "Actual Termination Factor: " << actualTerminationFactor << endl;
            cerr << endl;
            cerr << "=========================================" << endl;
            cerr << "Information on Evidence" << endl;
            cerr << "=========================================" << endl;
            cerr << "Skilling's log(Evidence): " << setprecision(12) << logEvidence << endl;
            cerr << "Keeton's log(Evidence): " << logMeanEvidence << endl;
            cerr << "Keeton's log(Live Evidence): " << logMeanLiveEvidence << endl;
            cerr << "Keeton's log(Total Evidence): " << logMeanTotalEvidence << endl;
            cerr << endl;
        }
        

        // If condition is verified, identify the clusters contained in the actual sample of nested objects

        if ((Niterations % NiterationsBeforeClustering) == 0)
        {
            
            Nclusters = clusterer.cluster(printOnTheScreen, nestedSampleOfParameters, clusterIndices);
        }


        // Evolve worst object with the new constraint logLikelihood > actualLogLikelihoodConstraint
        // Compute approximate sampling to find new point verifying the likelihood constraint

        ArrayXXd drawnSampleOfParameters = ArrayXXd::Zero(Ndimensions, 1);
        drawWithConstraint(nestedSampleOfParameters, Nclusters, clusterIndices, logTotalWidthInPriorMass, drawnSampleOfParameters, maxNdrawAttempts); 
       
        
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
        nestedSampleOfParameters.col(worst) = drawnSampleOfParameters;     // Update new set of parameters in the overall nested sample coordinates


        // Shrink prior mass interval
        
        logWidthInPriorMass -= 1.0 / Nobjects;

        
        // Update total width in prior mass from beginning to actual nested iteration

        logTotalWidthInPriorMass = Functions::logExpSum(logTotalWidthInPriorMass,logWidthInPriorMass);


        // Increase nested loop counter
        
        Niterations++;
    }
    while (terminationFactor < actualTerminationFactor);                          // Termination condition by Keeton 2011
    
 
    // Compute Skilling's error on the log(Evidence)
    
    logEvidenceError = sqrt(fabs(informationGain)/Nobjects);


    // Compute Keeton's errors on the log(Evidence)

    computeKeetonEvidenceError(printOnTheScreen, logMeanLiveEvidence);


    // Compute and print total computational time

    printComputationalTime(startTime);

}









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
}










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
}










// NestedSampler::getLogEvidenceError()
//
// PURPOSE:
//      Get private data member logEvidenceError.
//
// OUTPUT:
//      A double containing the Skilling's error on the logEvidence.
//

double NestedSampler::getLogEvidenceError()
{
    return logEvidenceError;
}










// NestedSampler::getLogMeanEvidence()
//
// PURPOSE:
//      Get private data member logMeanEvidence.
//
// OUTPUT:
//      A double containing the natural logarithm of the Keeton's mean evidence.
//

double NestedSampler::getLogMeanEvidence()
{
    return logMeanEvidence;
}









// NestedSampler::getLogMeanEvidenceError()
//
// PURPOSE:
//      Get private data member logMeanEvidenceError.
//
// OUTPUT:
//      A double containing the Keeton's error on the logMeanEvidence.
//

double NestedSampler::getLogMeanEvidenceError()
{
    return logMeanEvidenceError;
}










// NestedSampler::getLogMeanTotalEvidence()
//
// PURPOSE:
//      Get private data member logMeanTotalEvidence.
//
// OUTPUT:
//      A double containing the Keeton's mean total evidence.
//

double NestedSampler::getLogMeanTotalEvidence()
{
    return logMeanTotalEvidence;
}











// NestedSampler::getLogMeanTotalEvidenceError()
//
// PURPOSE:
//      Get private data member logMeanTotalEvidenceError.
//
// OUTPUT:
//      A double containing the Keeton's error on the logarithm of the mean total evidence.
//

double NestedSampler::getLogMeanTotalEvidenceError()
{
    return logMeanTotalEvidenceError;
}










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
}










// NestedSampler::getComputationalTime()
//
// PURPOSE:
//      Get private data member computationalTime.
//
// OUTPUT:
//      A double containing the final computational time of the process.
//

double NestedSampler::getComputationalTime()
{
    return computationalTime;
}











// NestedSampler::computeKeetonEvidenceError()
//
// PURPOSE:
//      Computes the total evidence error on the total evidence, according to the
//      description by Keeton C. 2011, MNRAS, 414,1418.
//
// INPUT:
//      printOnTheScreen: a boolean value specifying whether the results are to 
//      be printed on the screen or not.
//      logMeanLiveEvidence: a double containing the logarithm of the live evidence
//      compure by means of Keeton's equations.
//
// OUTPUT:
//      void
//

void NestedSampler::computeKeetonEvidenceError(const bool printOnTheScreen, const double logMeanLiveEvidence)
{
    double logAdditionalEvidenceError;
    double errorTerm1 = -DBL_MAX;
    double errorTerm2 = -DBL_MAX;
    double errorTerm3 = -DBL_MAX;
    double errorTerm4 = -DBL_MAX;
    double errorTerm5;

    for (int j = 0; j < Niterations; j++)
    {
        for (int i = 0; i < j; i++)
        {
            errorTerm2 = Functions::logExpSum(logLikelihoodOfPosteriorSample(i) + i*log(constant3), errorTerm2);
        }

        errorTerm1 = Functions::logExpSum(logLikelihoodOfPosteriorSample(j) + j*log(constant2) + errorTerm2, errorTerm1);
        errorTerm3 = Functions::logExpSum(logLikelihoodOfPosteriorSample(j) + j*log(constant2), errorTerm3);
        errorTerm2 = -DBL_MAX;
        errorTerm5 = Functions::logExpDifference(log(constant1) + j*log(constant3), j*log(constant2) - log(Nobjects));
        errorTerm4 = Functions::logExpSum(logLikelihoodOfPosteriorSample(j) + errorTerm5, errorTerm4);
    }
    
    logMeanEvidenceError = log(2.0) + log(constant1) - log(Nobjects) + errorTerm1;
    logMeanEvidenceError = Functions::logExpDifference(logMeanEvidenceError, -2.0*log(Nobjects) + 2.0*errorTerm3);
    
    logAdditionalEvidenceError = 2.0*logMeanLikelihoodOfLivePoints + Niterations*log(constant2) 
                                 + Functions::logExpDifference(Niterations*log(constant3), Niterations*log(constant3));
    logAdditionalEvidenceError = Functions::logExpSum(logAdditionalEvidenceError, log(2) + logMeanLikelihoodOfLivePoints + Niterations*log(constant2) + errorTerm4);
    logMeanTotalEvidenceError = sqrt(exp(logMeanEvidenceError) + exp(logAdditionalEvidenceError))/(exp(logMeanLiveEvidence) + exp(logMeanEvidence));
    logMeanEvidenceError = sqrt(exp(logMeanEvidenceError))/exp(logMeanEvidence);

    if (printOnTheScreen)
    {
        cerr << "=========================================" << endl;
        cout << "Final information on Evidence Uncertainty" << endl;
        cout << "=========================================" << endl;
        cerr << "Skilling's Uncertainty log(Evidence): " << logEvidenceError << endl; 
        cerr << "Keeton's Uncertainty log(Evidence): " << logMeanEvidenceError << endl; 
        cerr << "Keeton's Uncertainty log(Total Evidence): " << logMeanTotalEvidenceError << endl; 
        cerr << endl;
    }
}










// NestedSampler::printComputationalTime()
//
// PURPOSE:
//      Computes the total computational time of the nested sampling process
//      and prints the result expressed in either seconds, minutes or hours on the screen.
//
// INPUT:
//      startTime a double specifying the seconds at the moment the process started
//
// OUTPUT:
//      void
//

void NestedSampler::printComputationalTime(const double startTime)
{
    double endTime = time(0);
    computationalTime = endTime - startTime; 
    
    if (computationalTime < 60)
    {
        cerr << "=========================================" << endl;
        cerr << "Total Computational Time: " << computationalTime << " seconds" << endl;
        cerr << "=========================================" << endl;
    }
    else 
        if ((computationalTime >= 60) && (computationalTime < 60*60))
        {
            computationalTime = computationalTime/60.;
            cerr << "=========================================" << endl;
            cerr << "Total Computational Time: " << setprecision(3) << computationalTime << " minutes" << endl;
            cerr << "=========================================" << endl;
        }
    else 
        if (computationalTime >= 60*60)
        {
            computationalTime = computationalTime/(60.*60.);
            cerr << "=========================================" << endl;
            cerr << "Total Computational Time: " << setprecision(3) << computationalTime << " hours" << endl;
            cerr << "=========================================" << endl;
        }
}
