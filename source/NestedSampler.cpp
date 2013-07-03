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
  logEvidence(numeric_limits<double>::lowest())
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
//      terminationFactor:              a double specifying the amount of exceeding information to
//                                      terminate nested iterations. Default is 1.
//      NiterationsBeforeClustering:    number of nested iterations required to recompute
//                                      the clustering of the points in the sample.
//      maxNdrawAttempts:               the maximum number of attempts allowed when drawing from a single ellipsoid.
//
// OUTPUT:
//      void
//
// REMARKS: 
//      Eigen Matrices are defaulted column-major. Hence the 
//      nestedSample and posteriorSample are resized as 
//      (Ndimensions, ...), rather than (... , Ndimensions).
//

void NestedSampler::run(const double maxRatioOfRemainderToActualEvidence, const int NiterationsBeforeClustering, const int maxNdrawAttempts)
{
    int startTime = time(0);
    double logMeanLiveEvidence;
    double ratioOfRemainderToActualEvidence;


    // Set up the random number generator. It generates integers random numbers
    // between 0 and Nobjects-1, inclusive.

    uniform_int_distribution<int> discreteUniform(0, Nobjects-1);


    // Draw the initial sample from the prior PDF. Different coordinates of a point
    // can have different priors, so these have to be sampled individually.

    nestedSample.resize(Ndimensions, Nobjects);
    int beginIndex = 0;
    int NdimensionsOfCurrentPrior;
    ArrayXXd priorSample;

    for (int i = 0; i < ptrPriors.size(); i++)
    {
        // Some priors cover one particalar coordinate, others may cover two or more coordinates
        // Find out how many dimensions the current prior covers.

        NdimensionsOfCurrentPrior = ptrPriors[i]->getNdimensions();
        

        // Draw the subset of coordinates randomly from the current prior
        
        priorSample.resize(NdimensionsOfCurrentPrior, Nobjects);
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

    int Nclusters = 0;
    vector<int> clusterIndices(Nobjects);
    vector<int> clusterSizes;


    // Start the nested sampling loop. Each iteration, we'll replace the point with the worst likelihood.
    // New points are drawn from the prior, but with the constraint that they should have a likelihood
    // that is better than the currently worst one.

    do 
    {
        // Resize the arrays to make room for an additional point.
        // Do so without destroying the original contents.

        posteriorSample.conservativeResize(Ndimensions, Niterations + 1);  
        logLikelihoodOfPosteriorSample.conservativeResize(Niterations + 1);
        logWeightOfPosteriorSample.conservativeResize(Niterations + 1);
        

        // Find the point with the worst likelihood. This likelihood value will set a constraint
        // when drawing new points later on.
        
        int indexOfLivePointWithWorstLikelihood;
        worstLiveLogLikelihood = logLikelihood.minCoeff(&indexOfLivePointWithWorstLikelihood);
        double logWeight = logWidthInPriorMass + worstLiveLogLikelihood;                


        // Update the evidence Z and the information Gain
        
        double logEvidenceNew = Functions::logExpSum(logEvidence, logWeight);
        informationGain = exp(logWeight - logEvidenceNew) * worstLiveLogLikelihood 
                                        + exp(logEvidence - logEvidenceNew) * (informationGain + logEvidence) 
                                        - logEvidenceNew;
        logEvidence = logEvidenceNew;

        
        // Although we will replace the point with the worst likelihood in the live sample, we will save
        // it in our collection of posterior sample. Also save its likelihood value and its weight

        posteriorSample.col(Niterations) = nestedSample.col(indexOfLivePointWithWorstLikelihood); 
        logLikelihoodOfPosteriorSample(Niterations) = worstLiveLogLikelihood; 
        logWeightOfPosteriorSample(Niterations) = logWeight; 


        // Compute the (logarithm of) the mean likelihood of the set of live points.
        // Note that we are not computing mean(log(likelihood)) but log(mean(likelhood)).
        // Since we are only storing the log(likelihood) values, this results in a peculiar
        // way of computing the mean.
        
        logMeanLikelihoodOfLivePoints = logLikelihood(0);

        for (int m = 1; m < Nobjects; m++)
        {
            logMeanLikelihoodOfLivePoints = Functions::logExpSum(logMeanLikelihoodOfLivePoints, logLikelihood(m));
        }

        logMeanLikelihoodOfLivePoints -= log(Nobjects);
        
        
        // Compute mean live evidence

        logMeanLiveEvidence = logMeanLikelihoodOfLivePoints + Niterations * (log(Nobjects) - log(Nobjects + 1));


        // Compute the ratio of the evidence of the live sample to the actual Skilling's evidence.
        // Only when we gathered enough evidence, this ratio will be sufficiently small so that we can stop the iterations.

        ratioOfRemainderToActualEvidence = exp(logMeanLiveEvidence - logEvidence);
        

        // Find clusters in our live sample of points. Don't do this every iteration but only
        // every x iterations, where x is given by 'NiterationsBeforeClustering'.

        if ((Niterations % NiterationsBeforeClustering) == 0)
        {            
            Nclusters = clusterer.cluster(nestedSample, clusterIndices, clusterSizes, printOnTheScreen);
        }


        // Print current information on the screen, if required

        if (printOnTheScreen)
        {
            cerr << "Niter=" << Niterations << "  Nclusters:" << Nclusters 
            << "  WidthPriorMass= " << exp(logTotalWidthInPriorMass) 
            << "  EvidenceRatio=" << ratioOfRemainderToActualEvidence << "\r";
            cerr << "Skilling's log(Evidence): " << setprecision(12) << logEvidence << endl;
            cerr << "Information Gain: " << informationGain << endl;
        }



        // Draw a new point, which should replace the point with the worst likelihood.
        // This new point should be drawn from the prior, but with a likelihood greater 
        // than the current worst likelihood. The drawing algorithm may need a starting point,
        // for which we will take a randomly chosen point of the live sample (excluding the
        // worst point).

        int indexOfRandomlyChosenPoint = 0;
        if (Nobjects > 1)
        {
            // Select randomly an index of a sample point, but not the one of the worst point

            do 
            {
                // 0 <= indexOfRandomlyChosenPoint < Nobjects

                indexOfRandomlyChosenPoint = discreteUniform(engine);
            } 
            while (indexOfRandomlyChosenPoint == indexOfLivePointWithWorstLikelihood);
        }


        // drawnPoint will be a starting point as input, and will contain the newly drawn point as output

        ArrayXd drawnPoint = nestedSample.col(indexOfRandomlyChosenPoint); 
        double logLikelihoodOfDrawnPoint = 0.0;
        drawWithConstraint(nestedSample, Nclusters, clusterIndices, clusterSizes, 
                           logTotalWidthInPriorMass, drawnPoint, logLikelihoodOfDrawnPoint, maxNdrawAttempts); 


        // Replace the point with the worst likelihood with our newly drawn one.

        nestedSample.col(indexOfLivePointWithWorstLikelihood) = drawnPoint;
        logLikelihood(indexOfLivePointWithWorstLikelihood) = logLikelihoodOfDrawnPoint;
        

        // Shrink prior mass interval
        
        logWidthInPriorMass -= 1.0 / Nobjects;

        
        // Update total width in prior mass from beginning to actual nested iteration

        logTotalWidthInPriorMass = Functions::logExpSum(logTotalWidthInPriorMass, logWidthInPriorMass);


        // Increase nested loop counter
        
        Niterations++;
    }
    while (ratioOfRemainderToActualEvidence > maxRatioOfRemainderToActualEvidence);                          // Termination condition by Keeton 2011
    

    // If we get here, we sampled the parameter space well enough to gather enough evidence Z.
    // Add the remaining live sample of points to our collection of posterior points.

    int oldNpointsInPosterior = posteriorSample.cols();
    posteriorSample.conservativeResize(Ndimensions, oldNpointsInPosterior + Nobjects);          // First make enough room
    posteriorSample.block(0, oldNpointsInPosterior, Ndimensions, Nobjects) = nestedSample;      // Then copy the live sample to the posterior array
    

    // Compute Skilling's error on the log(Evidence)
    
    logEvidenceError = sqrt(fabs(informationGain)/Nobjects);


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
//      A double containing the natural logarithm of the Skilling's evidence.
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
    else 
        if (computationalTime >= 60*60*24)
        {
            computationalTime = computationalTime/(60.*60.*24.);
            cerr << "Total Computational Time: " << setprecision(3) << computationalTime << " days" << endl;
        }
}
