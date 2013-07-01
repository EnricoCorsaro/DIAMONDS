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
  uniform(0.0, 1.0),
  Niterations(0),
  informationGain(0.0), 
  logEvidence(-DBL_MAX),
  logMeanEvidence(-DBL_MAX),
  constant1(1.0/(Nobjects + 1)),
  constant2(constant1*Nobjects),
  constant3(1.0/((Nobjects + 2)*constant1))

{
   // Set the seed of the random generator using the clock

    clock_t clockticks = clock();
    engine.seed(clockticks);

    // The number of dimensions of the parameter space is the sum
    // of the dimensions covered by each of the priors

    Ndimensions = 0;
    
    for (int i = 0; i < ptrPriors.size(); i++)
    {
        // Get the number of dimensions from each type of prior

        Ndimensions += ptrPriors[i]->getNdimensions(); 
    }
} 









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
//      nestedSample and posteriorSample are resized as 
//      (Ndimensions, ...), rather than (... , Ndimensions).
//

void NestedSampler::run(const double maxRatioOfLiveToTotalEvidence, const int NiterationsBeforeClustering, const int maxNdrawAttempts)
{
    int startTime = time(0);
    double logMeanLiveEvidence;
    double ratioOfLiveToTotalEvidence;


    // Set up the random number generator. It generates integers random numbers
    // between 0 and Nobjects-1, inclusive.

    uniform_int_distribution<int> discreteUniform(0, Nobjects-1);


    // Draw the initial sample from the prior PDF. Different coordinates of a point
    // can have different priors, so these have to be sampled individually.

    nestedSample.resize(Ndimensions, Nobjects);
    int beginIndex = 0;

    for (int i = 0; i < ptrPriors.size(); i++)
    {
        // Some priors cover one particalar coordinate, others may cover two or more coordinates
        // Find out how many dimensions the current prior covers.

        const int NdimensionsOfCurrentPrior = ptrPriors[i]->getNdimensions();
        
        // Draw the subset of coordinates randomly from the current prior

        ArrayXXd priorSample(NdimensionsOfCurrentPrior, Nobjects);
        ptrPriors[i]->draw(priorSample);
        
        // Insert this random subset of coordinates into the total sample of coordinates of points

        nestedSample.block(beginIndex, 0, NdimensionsOfCurrentPrior, Nobjects) = priorSample;      
        
        // Move index to the beginning of the coordinate set of the next prior

        beginIndex += NdimensionsOfCurrentPrior;
    }



    // Compute the log(Likelihood) for each of our sample points

    logLikelihood.resize(Nobjects);
    
    for (int i = 0; i < Nobjects; ++i)
    {
        logLikelihood(i) = likelihood.logValue(nestedSample.col(i));
    }


    // Initialize the prior mass interval

    double logWidthInPriorMass = log(1.0 - exp(-1.0/Nobjects));
    double logTotalWidthInPriorMass = logWidthInPriorMass;
    
    
    // The nested sampling will involve finding clusters in the sample.
    // This will require the following containers.

    int Nclusters=0;
    vector<int> clusterIndices(Nobjects);
    vector<int> clusterSizes;


    // Start the nested sampling loop. Each iteration, we'll replace the point with the worse likelihood.
    // New points are drawn from the prior, but with the constraint that they should have a likelihood
    // that is better than the currently worse one.

    do 
    {
        // Resize the arrays to make room for an additional point.
        // Do so without destroying the original contents.

        posteriorSample.conservativeResize(Ndimensions, Niterations + 1);  
        logLikelihoodOfPosteriorSample.conservativeResize(Niterations + 1);
        logWeightOfPosteriorSample.conservativeResize(Niterations + 1);
        

        // Find the point with the worse likelihood. This likelihood value will set a constraint
        // when drawing new points later on.
        
        int indexOfLivePointWithWorseLikelihood;
        worseLiveLogLikelihood = logLikelihood.minCoeff(&indexOfLivePointWithWorseLikelihood);
        double logWeight = logWidthInPriorMass + worseLiveLogLikelihood;                


        // Update the evidence Z and the information Gain
        
        double logEvidenceNew = Functions::logExpSum(logEvidence, logWeight);
        informationGain = exp(logWeight - logEvidenceNew) * worseLiveLogLikelihood 
                                        + exp(logEvidence - logEvidenceNew) * (informationGain + logEvidence) 
                                        - logEvidenceNew;
        logEvidence = logEvidenceNew;

        
        // Although we will replace the point with the worse likelihood in the live sample, we will save
        // it in our collection of posterior sample. Also save its likelihood value and its weight

        posteriorSample.col(Niterations) = nestedSample.col(indexOfLivePointWithWorseLikelihood); 
        logLikelihoodOfPosteriorSample(Niterations) = worseLiveLogLikelihood; 
        logWeightOfPosteriorSample(Niterations) = logWeight; 


        // Compute the (logarithm of) the mean likelihood of the set of live points.
        // Note that we are not computing mean(log(likelihodd)) but log(mean(likelhood)).
        // Since we are only storing the log(likelihood) values, this results in a peculiar
        // way of computing the mean.
        
        logMeanLikelihoodOfLivePoints = logLikelihood(0);

        for (int m = 1; m < Nobjects; m++)
        {
            logMeanLikelihoodOfLivePoints = Functions::logExpSum(logMeanLikelihoodOfLivePoints, logLikelihood(m));
        }

        logMeanLikelihoodOfLivePoints -= log(Nobjects);
        
        
        // Compute mean evidence and mean live evidence

        logMeanLiveEvidence = logMeanLikelihoodOfLivePoints + Niterations * log(constant2);
        double logMeanEvidenceNew = Functions::logExpSum(worseLiveLogLikelihood + Niterations*log(constant2) - log(Nobjects), logMeanEvidence);
        logMeanEvidence = logMeanEvidenceNew;
        logMeanTotalEvidence = Functions::logExpSum(logMeanEvidence, logMeanLiveEvidence);


        // Compute the ratio of the evidence of the live sample to the evidence of the total sample (live + not-live).
        // Only when we gathered enough evidence, this ratio will be sufficiently small so that we can stop the iterations.

        ratioOfLiveToTotalEvidence = exp(logMeanLiveEvidence - logMeanEvidence);
        

        // Find clusters in our live sample of points. Don't do this every iteration but only
        // every x iterations, where x is given by 'NiterationsBeforeClustering'.

        if ((Niterations % NiterationsBeforeClustering) == 0)
        {            
            Nclusters = clusterer.cluster(nestedSample, clusterIndices, clusterSizes, printOnTheScreen);
        }


        // Print current information on the screen, if required

        if (printOnTheScreen)
        {
            cerr << "Niter=" << Niterations << "  Nclusters:" << Nclusters << "  WidthPriorMass= " << exp(logTotalWidthInPriorMass) << "  EvidenceRatio=" << ratioOfLiveToTotalEvidence << "\r";
            //cerr << "Skilling's log(Evidence): " << setprecision(12) << logEvidence << endl;
            //cerr << "Keeton's log(Evidence): " << logMeanEvidence << endl;
            //cerr << "Keeton's log(Live Evidence): " << logMeanLiveEvidence << endl;
            //cerr << "Keeton's log(Total Evidence): " << logMeanTotalEvidence << endl;
            //cerr << endl;
        }



        // Draw a new point, which should replace the point with the worse likelihood.
        // This new point should be drawn from the prior, but with a likelihood greater 
        // than the current worse likelihood. The drawing algorithm may need a starting point,
        // for which we will take a randomly chosen point of the live sample (excluding the
        // worse point).

        int indexOfRandomlyChosenPoint = 0;
        if (Nobjects > 1)
        {
            // Select randomly an index of a sample point, but not the one of the worse point

            do 
            {
                // 0 <= indexOfRandomlyChosenPoint < Nobjects

                indexOfRandomlyChosenPoint = discreteUniform(engine);
            } 
            while (indexOfRandomlyChosenPoint == indexOfLivePointWithWorseLikelihood);
        }


        // drawnPoint will be a starting point as input, and will contain the newly drawn point as output

        ArrayXd drawnPoint = nestedSample.col(indexOfRandomlyChosenPoint); 
        double logLikelihoodOfDrawnPoint = 0.0;
        drawWithConstraint(nestedSample, Nclusters, clusterIndices, clusterSizes, logTotalWidthInPriorMass, drawnPoint, logLikelihoodOfDrawnPoint, maxNdrawAttempts); 


        // Replace the point with the worse likelihood with our newly drawn one.

        nestedSample.col(indexOfLivePointWithWorseLikelihood) = drawnPoint;
        logLikelihood(indexOfLivePointWithWorseLikelihood) = logLikelihoodOfDrawnPoint;
        

        // Shrink prior mass interval
        
        logWidthInPriorMass -= 1.0 / Nobjects;

        
        // Update total width in prior mass from beginning to actual nested iteration

        logTotalWidthInPriorMass = Functions::logExpSum(logTotalWidthInPriorMass, logWidthInPriorMass);


        // Increase nested loop counter
        
        Niterations++;
    }
    while (ratioOfLiveToTotalEvidence > maxRatioOfLiveToTotalEvidence);                          // Termination condition by Keeton 2011
    
 
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
    
    cerr << endl;

    if (computationalTime < 60)
    {
        cerr << "Total Computational Time: " << computationalTime << " seconds" << endl;
    }
    else 
        if ((computationalTime >= 60) && (computationalTime < 60*60))
        {
            computationalTime = computationalTime/60.;
            cerr << "Total Computational Time: " << setprecision(3) << computationalTime << " minutes" << endl;
        }
    else 
        if (computationalTime >= 60*60)
        {
            computationalTime = computationalTime/(60.*60.);
            cerr << "Total Computational Time: " << setprecision(3) << computationalTime << " hours" << endl;
        }
}
