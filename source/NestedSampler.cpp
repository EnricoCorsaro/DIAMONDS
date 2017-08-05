#include "NestedSampler.h"


// NestedSampler::NestedSampler()
//
// PURPOSE: 
//      Constructor. Sets initial information, logEvidence and type 
//      of prior and likelihood distributions to be used. 
//
// INPUT:
//      printOnTheScreen:       Boolean value specifying whether the results are to 
//                              be printed on the screen or not.
//      initialNlivePoints:        Initial number of live points to start the nesting process
//      minNlivePoints:            Minimum number of live points allowed in the nesting process
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

NestedSampler::NestedSampler(const bool printOnTheScreen, const int initialNlivePoints, const int minNlivePoints, vector<Prior*> ptrPriors, 
                             Likelihood &likelihood, Metric &metric, Clusterer &clusterer)
: ptrPriors(ptrPriors),
  likelihood(likelihood),
  metric(metric),
  clusterer(clusterer),
  printOnTheScreen(printOnTheScreen),
  NlivePoints(initialNlivePoints),
  minNlivePoints(minNlivePoints),
  logCumulatedPriorMass(numeric_limits<double>::lowest()),
  logRemainingPriorMass(0.0),
  ratioOfRemainderToCurrentEvidence(numeric_limits<double>::max()),
  Niterations(0),
  updatedNlivePoints(initialNlivePoints),
  initialNlivePoints(initialNlivePoints),
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
//      logWeightOfPosteriorSample,logEvidenceOfPosteriorSample,logMeanLiveEvidenceOfPosteriorSample.
//
// INPUT:
//      livePointsReducer:                    An object of a class that takes care of the way the number of live points
//                                            is reduced within the nesting process
//      NinitialIterationsWithoutClustering:  The first N iterations, no clustering will happen. I.e. It will be assumed that
//                                            there is only 1 cluster containing all the points. This is often useful because 
//                                            initially the points may be sampled from a uniform prior, and we therefore don't 
//                                            expect any clustering before the algorithm is able to tune in on the island(s) of 
//                                            high likelihood. Clusters found in the first N initial iterations are therefore 
//                                            likely purely noise.
//      NiterationsWithSameClustering:        A new clustering will only happen every N iterations.
//      maxNdrawAttempts:                     The maximum number of attempts allowed when drawing from a single ellipsoid.
//      maxRatioOfRemainderToCurrentEvidence: The fraction of remainder evidence to gained evidence used to terminate 
//                                            the nested iteration loop. This value is also used as a tolerance on the final
//                                            evidence to update the number of live points in the nesting process.
//      pathPrefix:                           A string specifying the path where the output information from the Nested Sampler
//                                            has to be saved.
//
// OUTPUT:
//      void
//
// REMARKS: 
//      Eigen Matrices are defaulted column-major. Hence the nestedSample and posteriorSample are resized as 
//      (Ndimensions, ...), rather than (... , Ndimensions).
//

void NestedSampler::run(LivePointsReducer &livePointsReducer, const int NinitialIterationsWithoutClustering, 
                        const int NiterationsWithSameClustering, const int maxNdrawAttempts, 
                        const double maxRatioOfRemainderToCurrentEvidence, string pathPrefix)
{
    int startTime = time(0);
    double logMeanLiveEvidence;
    terminationFactor = maxRatioOfRemainderToCurrentEvidence;
    outputPathPrefix = pathPrefix;

    if (printOnTheScreen)
    {
        cerr << "------------------------------------------------" << endl;
        cerr << " Bayesian Inference problem has " << Ndimensions << " dimensions." << endl;
        cerr << "------------------------------------------------" << endl;
        cerr << endl;
    }


    // Save configuring parameters to an output ASCII file

    string fileName = "configuringParameters.txt";
    string fullPath = outputPathPrefix + fileName;
    File::openOutputFile(outputFile, fullPath);
   
   
    outputFile << "# List of configuring parameters used for the NSMC." << endl;
    outputFile << "# Row #1: Ndimensions" << endl;
    outputFile << "# Row #2: Initial(Maximum) NlivePoints" << endl;
    outputFile << "# Row #3: Minimum NlivePoints" << endl;
    outputFile << "# Row #4: NinitialIterationsWithoutClustering" << endl;
    outputFile << "# Row #5: NiterationsWithSameClustering" << endl;
    outputFile << "# Row #6: maxNdrawAttempts" << endl;
    outputFile << "# Row #7: terminationFactor" << endl;
    outputFile << "# Row #8: Niterations" << endl;
    outputFile << "# Row #9: Optimal Niterations" << endl;
    outputFile << "# Row #10: Final Nclusters" << endl;
    outputFile << "# Row #11: Final NlivePoints" << endl;
    outputFile << "# Row #12: Computational Time (seconds)" << endl;
    outputFile << Ndimensions << endl;
    outputFile << initialNlivePoints << endl;
    outputFile << minNlivePoints << endl;
    outputFile << NinitialIterationsWithoutClustering << endl;
    outputFile << NiterationsWithSameClustering << endl;
    outputFile << maxNdrawAttempts << endl;
    outputFile << terminationFactor << endl;


    // Set up the random number generator. It generates integer random numbers
    // between 0 and NlivePoints-1, inclusive.

    uniform_int_distribution<int> discreteUniform(0, NlivePoints-1);


    // Draw the initial sample from the prior PDF. Different coordinates of a point
    // can have different priors, so these have to be sampled individually.
    
    if (printOnTheScreen)
    {
        cerr << "------------------------------------------------" << endl;
        cerr << " Doing initial sampling of parameter space..." << endl;
        cerr << "------------------------------------------------" << endl;
        cerr << endl;
    }
        
    nestedSample.resize(Ndimensions, NlivePoints);
    int beginIndex = 0;
    int NdimensionsOfCurrentPrior;
    ArrayXXd priorSample;

    for (int i = 0; i < ptrPriors.size(); i++)
    {
        // Some priors cover one particalar coordinate, others may cover two or more coordinates
        // Find out how many dimensions the current prior covers.

        NdimensionsOfCurrentPrior = ptrPriors[i]->getNdimensions();
        

        // Draw the subset of coordinates randomly from the current prior
        
        priorSample.resize(NdimensionsOfCurrentPrior, NlivePoints);
        ptrPriors[i]->draw(priorSample);


        // Insert this random subset of coordinates into the total sample of coordinates of points

        nestedSample.block(beginIndex, 0, NdimensionsOfCurrentPrior, NlivePoints) = priorSample;      


        // Move index to the beginning of the coordinate set of the next prior

        beginIndex += NdimensionsOfCurrentPrior;
    }


    // Compute the log(Likelihood) for each of our points in the live sample

    logLikelihood.resize(NlivePoints);
    
    for (int i = 0; i < NlivePoints; ++i)
    {
        logLikelihood(i) = likelihood.logValue(nestedSample.col(i));
    }


    // Initialize the prior mass interval and cumulate it

    double logWidthInPriorMass = log(1.0 - exp(-1.0/NlivePoints));                                      // X_0 - X_1    First width in prior mass
    logCumulatedPriorMass = Functions::logExpSum(logCumulatedPriorMass, logWidthInPriorMass);           // 1 - X_1
    logRemainingPriorMass = Functions::logExpDifference(logRemainingPriorMass, logWidthInPriorMass);    // X_1


    // Initialize first part of width in prior mass for trapezoidal rule
    // X_0 = (2 - X_1), right-side boundary condition for trapezoidal rule

    double logRemainingPriorMassRightBound = Functions::logExpDifference(log(2), logRemainingPriorMass);    
    double logWidthInPriorMassRight = Functions::logExpDifference(logRemainingPriorMassRightBound,logRemainingPriorMass);


    // Find maximum log(Likelihood) value in the initial sample of live points. 
    // This information can be useful when reducing the number of live points adopted within the nesting process.

    logMaxLikelihoodOfLivePoints = logLikelihood.maxCoeff();


    // The nested sampling will involve finding clusters in the sample.
    // This will require the containers clusterIndices and clusterSizes.

    unsigned int Nclusters = 0;
    vector<int> clusterIndices(NlivePoints);        // clusterIndices must have the same number of elements as the number of live points
    vector<int> clusterSizes;                       // The number of live points counted in each cluster is updated everytime one live point
                                                    // is removed from the sample.


    // Start the nested sampling loop. Each iteration, we'll replace the point with the worst likelihood.
    // New points are drawn from the prior, but with the constraint that they should have a likelihood
    // that is better than the currently worst one.
    
    if (printOnTheScreen)
    {
        cerr << "-------------------------------" << endl;
        cerr << " Starting nested sampling...   " << endl;
        cerr << "-------------------------------" << endl;
        cerr << endl;
    }
        
    bool nestedSamplingShouldContinue = true;
    bool livePointsShouldBeReduced = (initialNlivePoints > minNlivePoints);       // Update live points only if required
    
    Niterations = 0;

    do 
    {
        // Resize the arrays to make room for an additional point.
        // Do so without destroying the original contents.

        posteriorSample.conservativeResize(Ndimensions, Niterations + 1);  
        logLikelihoodOfPosteriorSample.conservativeResize(Niterations + 1);
        logEvidenceOfPosteriorSample.conservativeResize(Niterations + 1);
        logMeanLiveEvidenceOfPosteriorSample.conservativeResize(Niterations + 1);
        logWeightOfPosteriorSample.conservativeResize(Niterations + 1);
        

        // Find the point with the worst likelihood. This likelihood value will set a constraint
        // when drawing new points later on.
        
        int indexOfLivePointWithWorstLikelihood;
        worstLiveLogLikelihood = logLikelihood.minCoeff(&indexOfLivePointWithWorstLikelihood);

        
        // Although we will replace the point with the worst likelihood in the live sample, we will save
        // it in our collection of posterior sample. Also save its likelihood value. The weight is 
        // computed and collected at the end of each iteration.

        posteriorSample.col(Niterations) = nestedSample.col(indexOfLivePointWithWorstLikelihood); 
        logLikelihoodOfPosteriorSample(Niterations) = worstLiveLogLikelihood; 


        // Compute the (logarithm of) the mean likelihood of the set of live points.
        // Note that we are not computing mean(log(likelihood)) but log(mean(likelhood)).
        // Since we are only storing the log(likelihood) values, this results in a peculiar
        // way of computing the mean. This will be used for computing the mean live evidence
        // at the end of the iteration.
        
        logMeanLikelihoodOfLivePoints = logLikelihood(0);

        for (int m = 1; m < NlivePoints; m++)
        {
            logMeanLikelihoodOfLivePoints = Functions::logExpSum(logMeanLikelihoodOfLivePoints, logLikelihood(m));
        }

        logMeanLikelihoodOfLivePoints -= log(NlivePoints);
                

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
                clusterSizes[0] = NlivePoints;
                fill(clusterIndices.begin(), clusterIndices.end(), 0);
            }
            else         
            {
                // After the first N initial iterations, we do a proper clustering.
                
                Nclusters = clusterer.cluster(nestedSample, clusterIndices, clusterSizes);
            }
        }


        // Draw a new point, which should replace the point with the worst likelihood.
        // This new point should be drawn from the prior, but with a likelihood greater 
        // than the current worst likelihood. The drawing algorithm may need a starting point,
        // for which we will take a randomly chosen point of the live sample (excluding the
        // worst point).

        int indexOfRandomlyChosenLivePoint = 0;
        
        if (NlivePoints > 1)
        {
            // Select randomly an index of a sample point, but not the one of the worst point

            do 
            {
                // 0 <= indexOfRandomlyChosenLivePoint < NlivePoints

                indexOfRandomlyChosenLivePoint = discreteUniform(engine);
            } 
            while (indexOfRandomlyChosenLivePoint == indexOfLivePointWithWorstLikelihood);
        }


        // drawnPoint will be a starting point as input, and will contain the newly drawn point as output

        ArrayXd drawnPoint = nestedSample.col(indexOfRandomlyChosenLivePoint);
        double logLikelihoodOfDrawnPoint = 0.0;
        bool newPointIsFound = drawWithConstraint(nestedSample, Nclusters, clusterIndices, clusterSizes, 
                                                  drawnPoint, logLikelihoodOfDrawnPoint, maxNdrawAttempts); 


        // If the adopted sampler produces an error (e.g. in the case of the ellipsoidal sampler a failure
        // in the ellipsoid matrix decomposition), then we can stop right here.
        
        nestedSamplingShouldContinue = verifySamplerStatus();
        if (!nestedSamplingShouldContinue) break;


        // If we didn't find a point with a better likelihood, then we can stop right here.
        
        if (!newPointIsFound)
        {
            nestedSamplingShouldContinue = false;
            cerr << "Can't find point with a better Likelihood." << endl; 
            cerr << "Stopping the nested sampling loop prematurely." << endl;
            break;
        }


        // Replace the point having the worst likelihood with our newly drawn one.

        nestedSample.col(indexOfLivePointWithWorstLikelihood) = drawnPoint;
        logLikelihood(indexOfLivePointWithWorstLikelihood) = logLikelihoodOfDrawnPoint;
       
        
        // If we got till here this is not the last iteration possible, hence 
        // update all the information for the next iteration. 
        // Check if the number of live points has not reached the minimum allowed,
        // and update it for the next iteration.

        if (livePointsShouldBeReduced)
        {
            // Update the number of live points for the current iteration based on the previous number.
            // If the number of live points reaches the minimum allowed 
            // then do not update the number anymore.

            updatedNlivePoints = livePointsReducer.updateNlivePoints();
            
            if (updatedNlivePoints > NlivePoints)
            {
                // Terminate program if new number of live points is greater than previous one
                    
                cerr << "Something went wrong in the reduction of the live points." << endl;
                cerr << "The new number of live points is greater than the previous one." << endl;
                cerr << "Quitting program. " << endl;
                break;
            }

                
            // If the lower bound for the number of live points has not been reached yet, 
            // the process should be repeated at the next iteration.
            // Otherwise the minimun number allowed is reached right now. In this case
            // stop the reduction process starting from the next iteration.
                
            livePointsShouldBeReduced = (updatedNlivePoints > minNlivePoints);

            if (updatedNlivePoints != NlivePoints)
            {
                // Resize all eigen arrays and vectors of dimensions NlivePoints according to 
                // new number of live points evaluated. In case previos and new number 
                // of live points coincide, no resizing is done.
                    
                vector<int> indicesOfLivePointsToRemove = livePointsReducer.findIndicesOfLivePointsToRemove(engine);

                    
                // At least one live point has to be removed, hence update the sample

                removeLivePointsFromSample(indicesOfLivePointsToRemove, clusterIndices, clusterSizes);
                        
                        
                // Since everything is fine update discreteUniform with the corresponding new upper bound

                uniform_int_distribution<int> discreteUniform2(0, updatedNlivePoints-1);
                discreteUniform = discreteUniform2;
            }
        }


        // Store the new number of live points in the vector containing this information.
        // This is done even if the new number is the same as the previous one.

        NlivePointsPerIteration.push_back(NlivePoints);

            
        // Compute the mean live evidence given the previous set of live points (see Keeton 2011, MNRAS) 

        logMeanLiveEvidence = logMeanLikelihoodOfLivePoints + Niterations * (log(NlivePoints) - log(NlivePoints + 1));


        // Compute the ratio of the evidence of the live sample to the current Skilling's evidence.
        // Only when we gathered enough evidence, this ratio will be sufficiently small so that we can stop the iterations.

        ratioOfRemainderToCurrentEvidence = exp(logMeanLiveEvidence - logEvidence);
        logEvidenceOfPosteriorSample(Niterations) = logEvidence;
        logMeanLiveEvidenceOfPosteriorSample(Niterations) = logMeanLiveEvidence;


        // Re-evaluate the stopping criterion, using the condition suggested by Keeton (2011)

        nestedSamplingShouldContinue = (ratioOfRemainderToCurrentEvidence > maxRatioOfRemainderToCurrentEvidence);


        // Shrink prior mass interval according to proper number of live points 
        // (see documentation by Enrico Corsaro October 2013). When reducing the number of live points 
        // the equation is a generalized version of that used by Skilling 2004. The equation
        // reduces to the standard case when the new number of live points is the same
        // as the previous one.

        // ---- Use the line below for simple rectangular rule ----
        // double logWeight = logWidthInPriorMass;
        // --------------------------------------------------------
        
        double logStretchingFactor = Niterations*((1.0/NlivePoints) - (1.0/updatedNlivePoints)); 
        logWidthInPriorMass = logRemainingPriorMass + Functions::logExpDifference(0.0, logStretchingFactor - 1.0/updatedNlivePoints);  // X_i - X_(i+1)

        
        // Compute the logWeight according to the trapezoidal rule 0.5*(X_(i-1) - X_(i+1)) 
        // and new contribution of evidence to be cumulated to the total evidence.
        // This is done in logarithmic scale by summing the right (X_(i-1) - X_i) and left part (X_i - X_(i+1)) 
        // of the total width in prior mass required for the trapezoidal rule. We do this computation at the end 
        // of the nested iteration because we need to know the new remaining prior mass of the next iteration.
            
        double logWidthInPriorMassLeft = logWidthInPriorMass; 

        
        // ---- Use the line below for trapezoidal rule ----

        double logWeight = log(0.5) + Functions::logExpSum(logWidthInPriorMassLeft, logWidthInPriorMassRight);
        double logEvidenceContributionNew = logWeight + worstLiveLogLikelihood;


        // Save log(Weight) of the current iteration

        logWeightOfPosteriorSample(Niterations) = logWeight;


        // Update the right part of the width in prior mass interval by replacing it with the left part

        logWidthInPriorMassRight = logWidthInPriorMass;


        // Print current information on the screen, if required

        if (printOnTheScreen)
        {
            if ((Niterations % 50) == 0)
            {
                cerr << "Nit: " << Niterations 
                     << "   Ncl: " << Nclusters 
                     << "   Nlive: " << NlivePoints
                     << "   CPM: " << exp(logCumulatedPriorMass)
                     << "   Ratio: " << ratioOfRemainderToCurrentEvidence
                     << "   log(E): " << logEvidence 
                     << "   IG: " << informationGain
                     << endl; 
            }
        }
        
        
        // Update the evidence and the information Gain
        
        double logEvidenceNew = Functions::logExpSum(logEvidence, logEvidenceContributionNew);
        informationGain = exp(logEvidenceContributionNew - logEvidenceNew) * worstLiveLogLikelihood 
                        + exp(logEvidence - logEvidenceNew) * (informationGain + logEvidence) 
                        - logEvidenceNew;
        logEvidence = logEvidenceNew;


        // Update total width in prior mass and remaining width in prior mass from beginning to current iteration
        // and use this information for the next iteration (if any)

        logCumulatedPriorMass = Functions::logExpSum(logCumulatedPriorMass, logWidthInPriorMass);
        logRemainingPriorMass = logStretchingFactor + logRemainingPriorMass - 1.0/updatedNlivePoints;


        // Update new number of live points in NestedSampler class 
            
        NlivePoints = updatedNlivePoints;


        // Increase nested loop counter
        
        Niterations++;
    }
    while (nestedSamplingShouldContinue);


    // Add the remaining live sample of points to our collection of posterior points 
    // (i.e parameter coordinates, likelihood values and weights)

    unsigned int oldNpointsInPosterior = posteriorSample.cols();

    posteriorSample.conservativeResize(Ndimensions, oldNpointsInPosterior + NlivePoints);          // First make enough room
    posteriorSample.block(0, oldNpointsInPosterior, Ndimensions, NlivePoints) = nestedSample;      // Then copy the live sample to the posterior array
    logWeightOfPosteriorSample.conservativeResize(oldNpointsInPosterior + NlivePoints);
    logWeightOfPosteriorSample.segment(oldNpointsInPosterior, NlivePoints).fill(logRemainingPriorMass - log(NlivePoints));  // Check if the best condition to impose 
    logLikelihoodOfPosteriorSample.conservativeResize(oldNpointsInPosterior + NlivePoints);
    logLikelihoodOfPosteriorSample.segment(oldNpointsInPosterior, NlivePoints) = logLikelihood; 


    // Compute Skilling's error on the log(Evidence)
    
    logEvidenceError = sqrt(fabs(informationGain)/NlivePoints);


    // Add Mean Live Evidence of the remaining live sample of points to the total log(Evidence) collected

    logEvidence = Functions::logExpSum(logMeanLiveEvidence, logEvidence);
    
    if (printOnTheScreen)
    {
        cerr << "------------------------------------------------" << endl;
        cerr << " Final log(E): " << logEvidence << " +/- " << logEvidenceError << endl;
        cerr << "------------------------------------------------" << endl;
    }

    // Print total computational time

    printComputationalTime(startTime);
    
    
    // Append information to existing output file and close stream afterwards
    
    outputFile << Niterations << endl;
    outputFile << static_cast<int>((NlivePoints*informationGain) + (NlivePoints*sqrt(Ndimensions*1.0))) << endl;
    outputFile << Nclusters << endl;
    outputFile << NlivePoints << endl;
    outputFile << computationalTime << endl;
}












// NestedSampler::removeLivePointsFromSample()
//
// PURPOSE:
//          Resizes all eigen arrays and vectors of dimensions NlivePoints according to 
//          new number of live points evaluated. The indices of the live points to be removed
//          are give as an input.
//          Also relative number of points in clusters are adjusted according to which live points
//          are removed.
//
// INPUT:   
//          indicesOfLivePointsToRemove:        A vector of integers containing the indices of the live points
//                                              that must be removed from the sample.
//          clusterIndices:                     A vector of integers containing the indices of the clusters
//                                              all the live points belong to
//          clusterSizes:                       A vector of integers containing the sizes of the clusters
// OUTPUT:
//      void
//

void NestedSampler::removeLivePointsFromSample(const vector<int> &indicesOfLivePointsToRemove, 
                                               vector<int> &clusterIndices, vector<int> &clusterSizes)
{
    int NlivePointsToRemove = indicesOfLivePointsToRemove.size();
    int NlivePointsAtCurrentIteration = clusterIndices.size();

    for (int m = 0; m < NlivePointsToRemove; ++m)
    {
        // Swap the last element of the set of live points with the chosen one 
        // and erase the last element. This is done for all the arrays that store information
        // about live points.
 
        ArrayXd nestedSamplePerLivePointCopy(Ndimensions);
        nestedSamplePerLivePointCopy = nestedSample.col(NlivePointsAtCurrentIteration-1);
        nestedSample.col(NlivePointsAtCurrentIteration-1) = nestedSample.col(indicesOfLivePointsToRemove[m]);
        nestedSample.col(indicesOfLivePointsToRemove[m]) = nestedSamplePerLivePointCopy;
        nestedSample.conservativeResize(Ndimensions, NlivePointsAtCurrentIteration-1);       
                
        double logLikelihoodCopy = logLikelihood(NlivePointsAtCurrentIteration-1);
        logLikelihood(NlivePointsAtCurrentIteration-1) = logLikelihood(indicesOfLivePointsToRemove[m]);
        logLikelihood(indicesOfLivePointsToRemove[m]) = logLikelihoodCopy;
        logLikelihood.conservativeResize(NlivePointsAtCurrentIteration-1);
        

        // In the case of clusterIndices also subtract selected live point from
        // corresponding clusterSizes in order to update the size of the cluster 
        // the live point belongs to.
                
        int clusterIndexCopy = clusterIndices[NlivePointsAtCurrentIteration-1];
        clusterIndices[NlivePointsAtCurrentIteration-1] = clusterIndices[indicesOfLivePointsToRemove[m]];
        --clusterSizes[clusterIndices[indicesOfLivePointsToRemove[m]]];
        clusterIndices[indicesOfLivePointsToRemove[m]] = clusterIndexCopy;
        clusterIndices.pop_back();

                
        // Reduce the current number of live points by one.
                
        --NlivePointsAtCurrentIteration;
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
   
    cerr << " Total Computational Time: ";

    if (computationalTime < 60)
    {
        cerr << computationalTime << " seconds" << endl;
    }
    else 
        if ((computationalTime >= 60) && (computationalTime < 60*60))
        {
            cerr << setprecision(3) << computationalTime/60. << " minutes" << endl;
        }
    else 
        if (computationalTime >= 60*60)
        {
            cerr << setprecision(3) << computationalTime/(60.*60.) << " hours" << endl;
        }
    else 
        if (computationalTime >= 60*60*24)
        {
            cerr << setprecision(3) << computationalTime/(60.*60.*24.) << " days" << endl;
        }
    
    cerr << "------------------------------------------------" << endl;
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

unsigned int NestedSampler::getNiterations()
{
    return Niterations;
}











// NestedSampler::getNdimensions()
//
// PURPOSE:
//      Get private data member Ndimensions.
//
// OUTPUT:
//      An integer containing the total number of
//      dimensions of the inference problem.
//

unsigned int NestedSampler::getNdimensions()
{
    return Ndimensions;
}












// NestedSampler::getNlivePoints()
//
// PURPOSE:
//      Get protected data member NlivePoints.
//
// OUTPUT:
//      An integer containing the current number of
//      live points.
//

int NestedSampler::getNlivePoints()
{
    return NlivePoints;
}












// NestedSampler::getInitialNlivePoints()
//
// PURPOSE:
//      Get protected data member initialNlivePoints.
//
// OUTPUT:
//      An integer containing the initial number of
//      live points.
//

int NestedSampler::getInitialNlivePoints()
{
    return initialNlivePoints;
}











// NestedSampler::getMinNlivePoints()
//
// PURPOSE:
//      Get protected data member minNlivePoints.
//
// OUTPUT:
//      An integer containing the minimum number of
//      live points allowed.
//

int NestedSampler::getMinNlivePoints()
{
    return minNlivePoints;
}












// NestedSampler::getLogCumulatedPriorMass()
//
// PURPOSE:
//      Get protected data member logCumulatedPriorMass.
//
// OUTPUT:
//      A double containing the natural logarithm of the cumulated prior mass.
//

double NestedSampler::getLogCumulatedPriorMass()
{
    return logCumulatedPriorMass;
}












// NestedSampler::getLogRemainingPriorMass()
//
// PURPOSE:
//      Get protected data member logRemainingPriorMass.
//
// OUTPUT:
//      A double containing the natural logarithm of the remaining prior mass.
//

double NestedSampler::getLogRemainingPriorMass()
{
    return logRemainingPriorMass;
}











// NestedSampler::getRatioOfRemainderToCurrentEvidence()
//
// PURPOSE:
//      Get protected data member ratioOfRemainderToCurrentEvidence.
//
// OUTPUT:
//      A double containing the ratio of the live evidence 
//      to the cumulated evidence.
//

double NestedSampler::getRatioOfRemainderToCurrentEvidence()
{
    return ratioOfRemainderToCurrentEvidence;
}












// NestedSampler::getLogMaxLikelihoodOfLivePoints()
//
// PURPOSE:
//      Get private data member logMaxLikelihoodOfLivePoints.
//
// OUTPUT:
//      A double containing the maximum log(Likelihood) value of the set of live points.
//

double NestedSampler::getLogMaxLikelihoodOfLivePoints()
{
    return logMaxLikelihoodOfLivePoints;
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











// NestedSampler::getTerminationFactor()
//
// PURPOSE:
//      Get private data member terminationFactor.
//
// OUTPUT:
//      A double containing the final value of the stopping condition for the nested process.
//

double NestedSampler::getTerminationFactor()
{
    return terminationFactor;
}











// NestedSampler::getNlivePointsPerIteration()
//
// PURPOSE:
//      Get protected data member NlivePointsPerIteration.
//
// OUTPUT:
//      A double containing the Skilling's error on the logEvidence.
//

vector<int> NestedSampler::getNlivePointsPerIteration()
{
    return NlivePointsPerIteration;
}













// NestedSampler::getNestedSample()
//
// PURPOSE:
//      Get private data member nestedSample.
//
// OUTPUT:
//      An eigen array containing the coordinates of the
//      current set of live points.
//

ArrayXXd NestedSampler::getNestedSample()
{
    return nestedSample;
}












// NestedSampler::getLogLikelihood()
//
// PURPOSE:
//      Get private data member logLikelihood.
//
// OUTPUT:
//      An eigen array containing the log(Likelihood) values of the
//      current set of live points.
//

ArrayXd NestedSampler::getLogLikelihood()
{
    return logLikelihood;
}











// NestedSampler::setLogEvidence()
//
// PURPOSE:
//      Set private data member logEvidence from the outside. 
//      Used when merging of the results is needed.
//
// OUTPUT:
//      void
//

void NestedSampler::setLogEvidence(double newLogEvidence)
{
    logEvidence = newLogEvidence;
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










// NestedSampler::setLogEvidenceError()
//
// PURPOSE:
//      Set private data member logEvidenceError from the outside. 
//      Used when merging of the results is needed.
//
// OUTPUT:
//      void
//

void NestedSampler::setLogEvidenceError(double newLogEvidenceError)
{
    logEvidenceError = newLogEvidenceError;
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










// NestedSampler::setInformationGain()
//
// PURPOSE:
//      Set private data member informationGain from the outside. 
//      Used when merging of the results is needed.
//
// OUTPUT:
//      void
//

void NestedSampler::setInformationGain(double newInformationGain)
{
    informationGain = newInformationGain;
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










// NestedSampler::setPosteriorSample()
//
// PURPOSE:
//      Set private data member posteriorSample from the outside. 
//      Used when merging of the results is needed.
//
// OUTPUT:
//      void
//

void NestedSampler::setPosteriorSample(ArrayXXd newPosteriorSample)
{
    Ndimensions = newPosteriorSample.rows();
    int Nsamples = newPosteriorSample.cols();
    posteriorSample.resize(Ndimensions, Nsamples);
    posteriorSample = newPosteriorSample;
}












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











// NestedSampler::setLogLikelihoodOfPosteriorSample()
//
// PURPOSE:
//      Set private data member logLikelihoodOfPosteriorSample from the outside. 
//      Used when merging of the results is needed.
//
// OUTPUT:
//      void
//

void NestedSampler::setLogLikelihoodOfPosteriorSample(ArrayXd newLogLikelihoodOfPosteriorSample)
{
    int Nsamples = newLogLikelihoodOfPosteriorSample.size();
    logLikelihoodOfPosteriorSample.resize(Nsamples);
    logLikelihoodOfPosteriorSample = newLogLikelihoodOfPosteriorSample;
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













// NestedSampler::setLogWeightOfPosteriorSample()
//
// PURPOSE:
//      Set private data member logWeightOfPosteriorSample from the outside. 
//      Used when merging of the results is needed.
//
// OUTPUT:
//      void
//

void NestedSampler::setLogWeightOfPosteriorSample(ArrayXd newLogWeightOfPosteriorSample)
{
    int Nsamples = newLogWeightOfPosteriorSample.size();
    logWeightOfPosteriorSample.resize(Nsamples);
    logWeightOfPosteriorSample = newLogWeightOfPosteriorSample;
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










// NestedSampler::getLogEvidenceOfPosteriorSample()
//
// PURPOSE:
//      Get private data member logEvidenceOfPosteriorSample.
//
// OUTPUT:
//      An eigen array containing the cumulated log(Evidence) values for each
//      nested iteration
//

ArrayXd NestedSampler::getLogEvidenceOfPosteriorSample()
{
    return logEvidenceOfPosteriorSample;
}











// NestedSampler::getLogMeanLiveEvidenceOfPosteriorSample()
//
// PURPOSE:
//      Get private data member logMeanLiveEvidenceOfPosteriorSample.
//
// OUTPUT:
//      An eigen array containing the log(MeanLiveEvidence) values remaining at each
//      nested iteration
//

ArrayXd NestedSampler::getLogMeanLiveEvidenceOfPosteriorSample()
{
    return logMeanLiveEvidenceOfPosteriorSample;
}










// NestedSampler::setOutputPathPrefix()
//
// PURPOSE:
//      Set private data member outputPathPrefix from the outside. 
//      Used when merging of the results is needed.
//
// OUTPUT:
//      void
//

void NestedSampler::setOutputPathPrefix(string newOutputPathPrefix)
{
    outputPathPrefix = newOutputPathPrefix;
}











// NestedSampler::getOutputPathPrefix()
//
// PURPOSE:
//      Get private data member outputPathPrefix.
//
// OUTPUT:
//      A string containing the full path of the output folder where all results
//      have to be saved.
//

string NestedSampler::getOutputPathPrefix()
{
    return outputPathPrefix;
}
