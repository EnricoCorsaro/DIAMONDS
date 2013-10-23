#include "NestedSampler.h"


// NestedSampler::NestedSampler()
//
// PURPOSE: 
//      Constructor. Sets initial information, logEvidence and type 
//      of prior and likelihood distributions to be used. 
//      Increases the number of live nested processes.
//
// INPUT:
//      printOnTheScreen:       Boolean value specifying whether the results are to 
//                              be printed on the screen or not.
//      initialNobjects:        Initial number of live points to start the nesting process
//      minNobjects:            Minimum number of live points allowed in the nesting process
//      ptrPriors:              Vector of pointers to Prior class objects
//      likelihood:             Likelihood class object used for likelihood sampling.
//      metric:                 Metric class object to contain the metric used in the problem.
//      clusterer:              Clusterer class object specifying the type of clustering algorithm to be used.
//
// REMARK:
//      The desired model for predictions is to be given initially to 
//      the likelihood object and is not feeded directly inside the 
//      nested sampling process.
//

NestedSampler::NestedSampler(const bool printOnTheScreen, const int initialNobjects, const int minNobjects, vector<Prior*> ptrPriors, 
                             Likelihood &likelihood, Metric &metric, Clusterer &clusterer)
: ptrPriors(ptrPriors),
  likelihood(likelihood),
  metric(metric),
  clusterer(clusterer),
  printOnTheScreen(printOnTheScreen),
  Nobjects(initialNobjects),
  logCumulatedPriorMass(numeric_limits<double>::lowest()),
  logRemainingPriorMass(0.0),
  minNobjects(minNobjects),
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
//      Destructor.
//

NestedSampler::~NestedSampler()
{
}









// NestedSampler::run()
//
// PURPOSE:
//      Start nested sampling computation. Save results in Eigen
//      Arrays logLikelihoodOfPosteriorSample, posteriorSample,
//      logWeightOfPosteriorSample.
//
// INPUT:
//      maxRatioOfRemainderToCurrentEvidence: The fraction of remainder evidence to gained evidence used to terminate 
//                                            the nested iteration loop. This value is also used as a tolerance on the final
//                                            evidence to update the number of live points in the nesting process.
//
//      NinitialIterationsWithoutClustering:  The first N iterations, no clustering will happen. I.e. It will be assumed that
//                                            there is only 1 cluster containing all the points. This is often useful because 
//                                            initially the points may be sampled from a uniform prior, and we therefore don't 
//                                            expect any clustering before the algorithm is able to tune in on the island(s) of 
//                                            high likelihood. Clusters found in the first N initial iterations are therefore 
//                                            likely purely noise.
//
//      NiterationsWithSameClustering:        A new clustering will only happen every N iterations.
//
//      maxNdrawAttempts:                     The maximum number of attempts allowed when drawing from a single ellipsoid.
//      livePointsReducer:                    An object of a class that takes care of the way the number of live points 
//                                            is reduced within the nesting process
//
// OUTPUT:
//      void
//
// REMARKS: 
//      Eigen Matrices are defaulted column-major. Hence the nestedSample and posteriorSample are resized as 
//      (Ndimensions, ...), rather than (... , Ndimensions).
//

void NestedSampler::run(LivePointsReducer &livePointsReducer, const double maxRatioOfRemainderToCurrentEvidence, const int NinitialIterationsWithoutClustering,
                        const int NiterationsWithSameClustering, const int maxNdrawAttempts)
{
    int startTime = time(0);
    double logMeanLiveEvidence;
    double ratioOfRemainderToCurrentEvidence;


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

    double logWidthInPriorMass = log(1.0 - exp(-1.0/Nobjects));         // First prior interval in the computation
    logCumulatedPriorMass = Functions::logExpSum(logCumulatedPriorMass,logWidthInPriorMass);
    
   
    // Evaluate max evidence contribution for the first iteration 

    logMaxLikelihoodOfLivePoints = logLikelihood.maxCoeff();
    logMaxEvidenceContribution = logMaxLikelihoodOfLivePoints;           // Initial remaining prior mass = 1


    // The nested sampling will involve finding clusters in the sample.
    // This will require the following containers.

    int Nclusters = 0;
    vector<int> clusterIndices(Nobjects);           // clusterIndices must have the same number of elements as the number of live points
    vector<int> clusterSizes;                       // The number of live points counted in each cluster is updated everytime one live point
                                                    // is removed from the sample.


    // Start the nested sampling loop. Each iteration, we'll replace the point with the worst likelihood.
    // New points are drawn from the prior, but with the constraint that they should have a likelihood
    // that is better than the currently worst one.

    bool nestedSamplingShouldContinue = true;
    bool livePointsShouldBeReduced = true;
    Niterations = 0;

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


        // Update the evidence and the information Gain
        
        double logEvidenceNew = Functions::logExpSum(logEvidence, logWeight);
        informationGain = exp(logWeight - logEvidenceNew) * worstLiveLogLikelihood 
                                        + exp(logEvidence - logEvidenceNew) * (informationGain + logEvidence) 
                                        - logEvidenceNew;
        logEvidence = logEvidenceNew;

        
        // Although we will replace the point with the worst likelihood in the live sample, we will save
        // it in our collection of posterior sample. Also save its likelihood value and its weight.

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
        
        
        // Compute mean live evidence (see Keeton 2011, MNRAS) 

        logMeanLiveEvidence = logMeanLikelihoodOfLivePoints + Niterations * (log(Nobjects) - log(Nobjects + 1));


        // Compute the ratio of the evidence of the live sample to the actual Skilling's evidence.
        // Only when we gathered enough evidence, this ratio will be sufficiently small so that we can stop the iterations.

        ratioOfRemainderToCurrentEvidence = exp(logMeanLiveEvidence - logEvidence);
        

        // Find clusters in our live sample of points. Don't do this every iteration but only
        // every x iterations, where x is given by 'NiterationsWithSameClustering'.
        
        if ((Niterations % NiterationsWithSameClustering) == 0)
        {            
            // Don't do clustering the first N iterations, where N is user-specified. That is, 
            // the first N iterations we assume that there is only 1 cluster containing all the points.
            // This is often useful because initially the points may be sampled from a uniform prior,
            // and we therefore don't expect any clustering _before_ the algorithm is able to tune in on 
            // the island(s) of high likelihood. Clusters found in the first N initial iterations are
            // therefore likely purely noise.
        
            if (Niterations < NinitialIterationsWithoutClustering)
            {
                // There is only 1 cluster, containing all objects. All points have the same cluster
                // index, namely 0.
                       
                Nclusters = 1;
                clusterSizes.resize(1);
                clusterSizes[0] = Nobjects;
                fill(clusterIndices.begin(), clusterIndices.end(), 0);
            }
            else         
            {
                // After the first N initial iterations, we do a proper clustering.
                
                Nclusters = clusterer.cluster(nestedSample, clusterIndices, clusterSizes, printOnTheScreen);
            }
        }


        // Print current information on the screen, if required

        if (printOnTheScreen)
        {
            if ((Niterations % 50) == 0)
            {
                cerr << "Nit: " << Niterations << "   Ncl: " << Nclusters 
                     << "   Nlive: " << Nobjects
                     << "   CPM: " << exp(logCumulatedPriorMass)
                     << "   Ratio: " << ratioOfRemainderToCurrentEvidence
                     << "   log(E): " << logEvidence 
                     << "   IG: " << informationGain
                     << endl;
            }
        }


        // Draw a new point, which should replace the point with the worst likelihood.
        // This new point should be drawn from the prior, but with a likelihood greater 
        // than the current worst likelihood. The drawing algorithm may need a starting point,
        // for which we will take a randomly chosen point of the live sample (excluding the
        // worst point).

        int indexOfRandomlyChosenLivePoint = 0;
        
        if (Nobjects > 1)
        {
            // Select randomly an index of a sample point, but not the one of the worst point

            do 
            {
                // 0 <= indexOfRandomlyChosenLivePoint < Nobjects

                indexOfRandomlyChosenLivePoint = discreteUniform(engine);
            } 
            while (indexOfRandomlyChosenLivePoint == indexOfLivePointWithWorstLikelihood);
        }


        // drawnPoint will be a starting point as input, and will contain the newly drawn point as output

        ArrayXd drawnPoint = nestedSample.col(indexOfRandomlyChosenLivePoint);
        double logLikelihoodOfDrawnPoint = 0.0;
        bool newPointIsFound = drawWithConstraint(nestedSample, Nclusters, clusterIndices, clusterSizes, 
                                                  drawnPoint, logLikelihoodOfDrawnPoint, maxNdrawAttempts); 

        // If we didn't find a point with a better likelihood, then we can stop right here.
        
        if (!newPointIsFound)
        {
            nestedSamplingShouldContinue = false;
            cerr << "Can't find point with a better Likelihood" << endl; 
            cerr << "Stopping the nested sampling loop prematurely." << endl;
            break;
        }


        // Replace the point having the worst likelihood with our newly drawn one.

        nestedSample.col(indexOfLivePointWithWorstLikelihood) = drawnPoint;
        logLikelihood(indexOfLivePointWithWorstLikelihood) = logLikelihoodOfDrawnPoint;
        
        
        // Increase nested loop counter
        
        Niterations++;


        // Re-evaluate the stopping criterion, using the condition suggested by Keeton (2011)

        nestedSamplingShouldContinue = (ratioOfRemainderToCurrentEvidence > maxRatioOfRemainderToCurrentEvidence);
       
        
        // If the current iteration is not the last iteratiom and 
        // the number of live points has not reached the minimum allowed,
        // then update the number of live points for the next iteration.

        if (nestedSamplingShouldContinue && livePointsShouldBeReduced)
        {
            // Evaluate max evidence contribution for the next iteration 

            logMaxLikelihoodOfLivePoints = logLikelihood.maxCoeff();
            double logMaxEvidenceContributionNew = logMaxLikelihoodOfLivePoints + logRemainingPriorMass;


            // Number of live points for the actual iteration

            int NobjectsAtCurrentIteration = Nobjects;


            // Update the number of live points for the next iteration. If the number of live points reaches the minimum allowed
            // then do not update their number anymore.

            livePointsShouldBeReduced = updateNobjects(logMaxEvidenceContributionNew, maxRatioOfRemainderToCurrentEvidence);
    
            
            // Number of live points to be removed

            int NobjectsToRemove = NobjectsAtCurrentIteration - Nobjects;                
           

            // Resize all eigen arrays of dimensions Nobjects according to new number of live points evaluated.
            // Pick one live point randomly and Remove it. Repeat the process up to the total number of live points to be removed.

            ArrayXd nestedSamplePerLivePointCopy(Ndimensions);

            for (int m = 0; m < NobjectsToRemove; ++m)
            {
                // Rescale uniform random integer generator with new number of live points

                uniform_int_distribution<int> discreteUniform2(0, NobjectsAtCurrentIteration-1);
                

                // Select randomly one live point from the actual sample

                int indexOfLivePointToRemove = discreteUniform2(engine);


                // Swap the last element of the set of live points with the chosen one 
                // and erase the last element. This is done for all the arrays that store information
                // about live points.

                nestedSamplePerLivePointCopy = nestedSample.col(NobjectsAtCurrentIteration-1);
                nestedSample.col(NobjectsAtCurrentIteration-1) = nestedSample.col(indexOfLivePointToRemove);
                nestedSample.col(indexOfLivePointToRemove) = nestedSamplePerLivePointCopy;
                nestedSample.conservativeResize(Ndimensions, NobjectsAtCurrentIteration-1);       
                
                double logLikelihoodCopy = logLikelihood(NobjectsAtCurrentIteration-1);
                logLikelihood(NobjectsAtCurrentIteration-1) = logLikelihood(indexOfLivePointToRemove);
                logLikelihood(indexOfLivePointToRemove) = logLikelihoodCopy;
                logLikelihood.conservativeResize(NobjectsAtCurrentIteration-1);
              

                // In the case of clusterIndices also subtract selected live point from 
                // corresponding clusterSizes in order to update the size of the cluster to which
                // the live point belongs.
                
                int clusterIndexCopy = clusterIndices[NobjectsAtCurrentIteration-1];
                clusterIndices[NobjectsAtCurrentIteration-1] = clusterIndices[indexOfLivePointToRemove];
                --clusterSizes[clusterIndices[indexOfLivePointToRemove]];
                clusterIndices[indexOfLivePointToRemove] = clusterIndexCopy;
                clusterIndices.pop_back();

                
                // Reduce the actual number of live points by one.
                
                --NobjectsAtCurrentIteration;
            }


            // Update discreteUniform with final number of live points

            uniform_int_distribution<int> discreteUniform3(0, Nobjects-1);
            discreteUniform = discreteUniform3;
        }
        

        // Shrink prior mass interval according to proper number of live points
        
        logWidthInPriorMass -= 1.0 / Nobjects;

        
        // Update total width in prior mass and remaining width in prior mass from beginning to actual iteration
        // and use this information for the next iteration (if any)

        logCumulatedPriorMass = Functions::logExpSum(logCumulatedPriorMass, logWidthInPriorMass);
        logRemainingPriorMass = log(1.0 - exp(logCumulatedPriorMass));
    }
    while (nestedSamplingShouldContinue);  

    
    // Add the remaining live sample of points to our collection of posterior points 
    // (i.e parameter coordinates, likelihood values and weights)

    int oldNpointsInPosterior = posteriorSample.cols();
    
    posteriorSample.conservativeResize(Ndimensions, oldNpointsInPosterior + Nobjects);          // First make enough room
    posteriorSample.block(0, oldNpointsInPosterior, Ndimensions, Nobjects) = nestedSample;      // Then copy the live sample to the posterior array
    logWeightOfPosteriorSample.conservativeResize(oldNpointsInPosterior + Nobjects);
    logWeightOfPosteriorSample.segment(oldNpointsInPosterior, Nobjects) = logWidthInPriorMass + logLikelihood;                
    logLikelihoodOfPosteriorSample.conservativeResize(oldNpointsInPosterior + Nobjects);
    logLikelihoodOfPosteriorSample.segment(oldNpointsInPosterior, Nobjects) = logLikelihood; 


    // Compute Skilling's error on the log(Evidence)
    
    logEvidenceError = sqrt(fabs(informationGain)/Nobjects);


    // Add Mean Live Evidence of the remaining live sample of points to the total log(Evidence) collected

    logEvidence = Functions::logExpSum(logMeanLiveEvidence, logEvidence);


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













// NestedSampler::updateNobjects()
//
// PURPOSE:
//      Computes the updated number of live points based on the gaining of evidence
//      achieved at the given iteration. The equation adopted was introduced
//      by Feroz F. et al. (2009), MNRAS, 398, 1601.
//
// INPUT:
//      logMaxEvidenceContributionNew:          The updated maximum contribution of the evidence
//                                              expressed in natural logarithm.
//      maxRatioOfRemainderToCurrentEvidence:    The fraction of remainder evidence to gained evidence used to terminate 
//                                              the nested iteration loop. This value is also used as a tolerance on the final
//                                              evidence to update the number of live points in the nesting process.
//
// OUTPUT:
//      An integer containing the new number of live points to be used.
//

bool NestedSampler::updateNobjects(double logMaxEvidenceContributionNew, double maxRatioOfRemainderToCurrentEvidence)
{
    // Evaluate the new number of live points to be used in the next iteration of the nesting loop
    
    int updatedNobjects = static_cast<int>(Nobjects - minNobjects * 
                          exp(Functions::logExpDifference(logMaxEvidenceContribution,logMaxEvidenceContributionNew)) /
                          (exp(logMaxEvidenceContributionNew)*(1 - maxRatioOfRemainderToCurrentEvidence)));
    
    assert(updatedNobjects <= Nobjects);

    if (updatedNobjects > minNobjects)
    {
        // If minimum number of live points allowed has not been reached, 
        // save the new number and use it for the next nested iteration

        Nobjects = updatedNobjects;
        return true;
    }
    else











// NestedSampler::getPosteriorSample()
//
// PURPOSE:
//      Get private data member posteriorSample.
//
// OUTPUT:
//      An eigen array containing the coordinates of the
//      final posterior sample.
//

ArrayXXd NestedSampler::getPosteriorSample()
{
    return posteriorSample;
}












// NestedSampler::getLogLikelihoodOfPosteriorSample()
//
// PURPOSE:
//      Get private data member logLikelihoodOfPosteriorSample.
//
// OUTPUT:
//      An eigen array containing the log(Likelihood) values of the
//      final posterior sample.
//

ArrayXd NestedSampler::getLogLikelihoodOfPosteriorSample()
{
    return logLikelihoodOfPosteriorSample;
}














// NestedSampler::getLogWeightOfPosteriorSample()
//
// PURPOSE:
//      Get private data member logWeightOfPosteriorSample.
//
// OUTPUT:
//      An eigen array containing the log(Weight) values of the
//      final posterior sample.
//

ArrayXd NestedSampler::getLogWeightOfPosteriorSample()
{
    return logWeightOfPosteriorSample;
}















    {
        // Otherwise continue the nesting process by using 
        // the minimum number of live points allowed

        return false;
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

