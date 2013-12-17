#include "Results.h"


// Results::Results()
//
// PURPOSE: 
//      Constructor. Sets nested sampler object and output directory path.
//
// INPUT:
//      nestedSampler:  a NestedSampler class object used as the container of
//                      information to write from.
//

Results::Results(NestedSampler &nestedSampler)
: nestedSampler(nestedSampler)
{
} 










// Results::~Results()
//
// PURPOSE: 
//      Destructor.
//

Results::~Results()
{

}












// Results::posteriorProbability()
//
// PURPOSE:
//      Saves the posterior probability for the sample, obtained by
//      applying the Bayes theorem to the information coming from the 
//      nested sampling, into a one-dimensional Eigen Array. 
//
// OUTPUT:
//      An Eigen Array containing the values of the posterior probability
//      sorted according to the nesting process. 
// 
// REMARK:
//      Values are probabilities (and not probability densities),
//      hence their sum must equal the unity.
//      The array computed according to the definition of posterior probability
//      from the nesting algorithm is finally normalized by the sum
//      of its elements. This allows to get rid of small deviations
//      caused by the approximated value of the evidence.
//

ArrayXd Results::posteriorProbability()
{
    // Apply Bayes Theorem in logarithmic expression

    ArrayXd logPosteriorDistribution = nestedSampler.getLogWeightOfPosteriorSample() + 
                                       nestedSampler.getLogLikelihoodOfPosteriorSample() - 
                                       nestedSampler.getLogEvidence();
    ArrayXd posteriorDistribution = logPosteriorDistribution.exp();
   
   
    // Since evidence is approximate, ensure that the sum of all probabilities equals the unity
    // when returning the array.

    return posteriorDistribution/posteriorDistribution.sum();
}











// Results:writeMarginalDistributionToFile()
//
// PURPOSE:
//      Writes the interpolated grid of parameter values and its 
//      corresponding marginal distribution as derived from parameterEstimation()
//      in an output ASCII file of two-columns format. 
//      This information can then be used for plotting the marginal distribution  
//      and for computing further operations on it, such as the derivation of error bars on the free parameter.
//
// INPUT:
//      pathPrefix:                 a string variable containing the desired path for saving the output file.
//      parameterNumber:            an integer containing the number of the free parameter to save, together with
//                                  its marginal distribution.
//      
// OUTPUT:
//      void.
// 

void Results::writeMarginalDistributionToFile(string pathPrefix, const int parameterNumber)
{
    // Input arrays must contain the same number of elements 

    int Nrows = parameterValuesInterpolated.size();
    assert(Nrows == marginalDistributionInterpolated.size());

    ArrayXXd parameterDistribution(Nrows,2);
    parameterDistribution.col(0) = parameterValuesInterpolated;
    parameterDistribution.col(1) = marginalDistributionInterpolated;


    // Include the row number with preceding zeros in the filename
        
    ostringstream numberString;
    numberString << setfill('0') << setw(3) << parameterNumber;
    string outputFileName = "marginal";
    string fullPath = pathPrefix + outputFileName + numberString.str() + ".txt";
    
    
    // Write the two-columns array in an ASCII file

    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
    outputFile << "# Marginal distribution of Akima-interpolated points from nested sampling." << endl;
    outputFile << "# Column #1: Parameter values" << endl;
    outputFile << "# Column #2: Marginal distribution values (probability only)" << endl;
    outputFile << scientific << setprecision(9);
    File::arrayXXdToFile(outputFile, parameterDistribution);
    outputFile.close();
}












// Results:computeCredibleLimits()
//
// PURPOSE:
//      Computes the shortest Bayesian credible intervals (CI) from a marginal distribution.
//
// INPUT:
//      credibleLevel:              a double number providing the desired credible 
//                                  level to be computed.
//      Nbins:                      an integer specifying the number of points of the distribution to interpolate
//      NinterpolationsPerBin:      an integer containing the number of desired points to interpolate between.
//                                  two consecutive input data points.
//      
// OUTPUT:
//      A one-dimensional Eigen Array containing the lower and upper credible limits, i.e. the boundaries
//      of the shortest credible intervals identified.
//
// REMARK:
//      The input number of interpolated points is to be considered as an estimate of the real number of interpolations occurring.
//      This is because the number of interpolations also depends on the spacing of the input grid for each bin, which is
//      not required to be regular.
//

ArrayXd Results::computeCredibleLimits(const double credibleLevel, const int Nbins, const int NinterpolationsPerBin)
{
    // Interpolate rebinned marginal distribution by using a NinterpolationsPerBin times finer grid. 
    // This allows for a better computation of the credible intervals. Note that while the rebinned 
    // marginal distribution has still sum equal to 1, the interpolated marginal distribution must 
    // be renormalized after interpolation. This will scale down the amplitudes of a factor of the order of NinterpolationsPerBin. 
    // Such a rescaling has no effect on the computation of the credible intervals,
    // which instead depend upon the relative variation from point to point in the distribution.

    int Ninterpolations = Nbins*NinterpolationsPerBin;


    // Take the average bin width of the rebinned data for computing the bin width of the interpolated data. This information is more
    // realistic as it better resembles the original sampling from the nesting process.

    double binWidth = (parameterValuesRebinned.segment(1, Nbins-1) - parameterValuesRebinned.segment(0, Nbins-1)).sum() / (Nbins*1.0);
    double interpolatedBinWidth = binWidth/(NinterpolationsPerBin*1.0);
    double parameterMinimumRebinned = parameterValuesRebinned.minCoeff();

    
    // Initialize the abscissa of the interpolated distribution
    
    parameterValuesInterpolated.resize(Ninterpolations);

    for (int j = 0; j < Ninterpolations; ++j)
    {
        parameterValuesInterpolated(j) = parameterMinimumRebinned + j*interpolatedBinWidth;
    }


    // Apply a cubic spline interpolation for the new interpolated values of the marginal distribution

    marginalDistributionInterpolated = Functions::cubicSplineInterpolation(parameterValuesRebinned, marginalDistributionRebinned, parameterValuesInterpolated);


    // Check if negative values (even if small) are present and if so, set them to zero.
    // This may occur in the case of small fluctuations around zero from the cubic spline. Probabilities cannot be negative.

    vector<int> selectedIndices = Functions::findArrayIndicesWithinBoundaries(marginalDistributionInterpolated, -1e99, -1e-99);
    
    for (int i = 0; i < selectedIndices.size(); ++i)
    {
        marginalDistributionInterpolated(selectedIndices[i]) = 0.0;
    }
        

    // Renormalize interpolated marginal distribution by the sum of its elements

    marginalDistributionInterpolated /= marginalDistributionInterpolated.sum();


    // Resize corresponding abscissa because total number of interpolated bins might be changed after interpolation (e.g. if truncated)

    if (marginalDistributionInterpolated.size() != Ninterpolations)
    {
        Ninterpolations = marginalDistributionInterpolated.size();
        parameterValuesInterpolated.conservativeResize(Ninterpolations);
    }


    // Compute the "shortest" credible intervals (CI) and save the corresponding credible limits

    int max = 0;
    marginalDistributionMode = marginalDistributionInterpolated.maxCoeff(&max);
    ArrayXd marginalDistributionLeft(max + 1);       // Marginal distribution up to mode value (included)
    ArrayXd parameterValuesLeft(max + 1);            // Parameter range up to mode value (included)
    double limitProbabilityRight = marginalDistributionInterpolated(max);
    double limitParameterRight = parameterValuesInterpolated(max);
    double limitProbabilityLeft = limitProbabilityRight;
    double limitParameterLeft = limitParameterRight;


    // Difference distribution to find point belonging to the distribution to the left
    // of the mode having closest probability to that identified in the part of the
    // distribution to the right of the mode.

    ArrayXd marginalDifferenceLeft = ArrayXd::Zero(max + 1);          
                                                

    // Copy left-hand part of total marginal distributon into separate array.
    // Also copy the corresponding parameter values inro a separate array.
    // Since parameterValues is already sorted in ascending order this can be done quickly.

    marginalDistributionLeft = marginalDistributionInterpolated.segment(0, max + 1);
    parameterValuesLeft = parameterValuesInterpolated.segment(0, max + 1);


    // Count number of steps (bins) in the distribution to the right side of its modal value.
    // Define the limit probability on the right-hand distribution 
    // and the corresponding value of the parameter.
    // Then compute the difference of the left-hand distribution with the right limit probability.
    // Finally find the point in the left part having closest probability to the right limit,
    // save its probability and the corresponding value of the parameter.
    // Cumulate the probability within the range and repeat the process until 
    // probability >= probability(credibleLevel).

    int stepRight = 0;          
    int min = 0;
    double credibleLevelFraction = credibleLevel/100.;
    double totalProbability = 0.0;

    while (totalProbability < (credibleLevelFraction) && (max + stepRight < Ninterpolations))             
    {
        totalProbability = 0.0;


        // Find the probability and the corresponding value of the parameter at the right edge of the interval

        limitProbabilityRight = marginalDistributionInterpolated(max + stepRight);  
        limitParameterRight = parameterValuesInterpolated(max + stepRight);         

            
        // Find which point in the left part is the closest in probability to that of the right edge.

        marginalDifferenceLeft = (marginalDistributionLeft - limitProbabilityRight).abs();                                      
        limitProbabilityLeft = marginalDifferenceLeft.minCoeff(&min);       


        // Save the point and its probability as the left edge of the search interval

        limitProbabilityLeft = marginalDistributionLeft(min);           
        limitParameterLeft = parameterValuesLeft(min);                 


        // Cumulate the probability within the identified interval

        int intervalSize = max + stepRight + 1 - min;
        totalProbability = marginalDistributionInterpolated.segment(min, intervalSize).sum();


        // Move one step further to the right

        ++stepRight;
    }

    ArrayXd credibleLimits(2);
    credibleLimits << limitParameterLeft, limitParameterRight; 
    
    return credibleLimits;
}














// Results:parameterEstimation()
//
// PURPOSE:
//      Computes the expectation, median and mode values from the 
//      marginalized posterior probability. 
//      Shortest Bayesian credible intervals (CI) are also computed.
//
// INPUT:
//      credibleLevel:  a double number providing the desired credible 
//                      level to be computed. Default value corresponds to
//                      credible level of 68.27 %.
//      
// OUTPUT:
//      A bidimensional Eigen Array containing all the estimators of the
//      free parameters (one parameter per row). The output order is:
//      (1) Mean
//      (2) Median
//      (3) Mode
//      (4) II Moment
//      (5) Lower CL
//      (6) Upper CL
// 

ArrayXXd Results::parameterEstimation(const double credibleLevel)
{
    int Ndimensions = nestedSampler.getPosteriorSample().rows();
    ArrayXd posteriorDistribution = posteriorProbability();
    assert (nestedSampler.getPosteriorSample().cols() == posteriorDistribution.size());
    ArrayXXd parameterEstimates(Ndimensions, 6);


    // Loop over all free parameters

    for (int i = 0; i < Ndimensions; ++i)
    {
        ArrayXd parameterValues = nestedSampler.getPosteriorSample().row(i);
        ArrayXd marginalDistribution = posteriorDistribution;


        // Sort elements of array parameterValues in ascending
        // order and sort elements of array marginalDistribution accordingly.
        // Use a mergesort algorithm for more efficiency with a large number of array elements.
        
        Functions::topDownMergeSort(parameterValues, marginalDistribution); 
        

        // Compute the mean value (expectation value, or first moment) of the sample distribution
       
        double parameterMean = (parameterValues * marginalDistribution).sum();
        parameterEstimates(i,0) = parameterMean;


        // Compute second moment of the sample distribution (i.e. the standard deviation in case of a normal distribution)

        double secondMoment = ((parameterValues - parameterMean).pow(2) * marginalDistribution).sum();

        
        // Compute the median value (value corresponding to 50% of probability)

        int k = 0;
        double totalProbability = 0.0;
        double parameterMedian = parameterValues(0);
        
        while (totalProbability < 0.50)
        {
            parameterMedian = parameterValues(k);
            totalProbability += marginalDistribution(k);
            k++;
        }
        
        parameterEstimates(i,1) = parameterMedian;


        // Find the mode value (parameter corresponding to maximum probability value)

        int max = 0;                                    // Subscript corresponding to mode value
        double marginalDistributionMode = marginalDistribution.maxCoeff(&max);
        double parameterMode = parameterValues(max);
        parameterEstimates(i,2) = parameterMode;

       
        // Save the second moment of the distribution
       
        parameterEstimates(i,3) = secondMoment;


        // Compute the "shortest" credible intervals (CI) and save the corresponding credible limits

        marginalDistributionMode = marginalDistribution.maxCoeff(&max);
        ArrayXd marginalDistributionLeft(max + 1);       // Marginal distribution up to mode value (included)
        ArrayXd parameterValuesLeft(max + 1);            // Parameter range up to mode value (included)
        double limitProbabilityRight = marginalDistribution(max);
        double limitParameterRight = parameterValues(max);
        double limitProbabilityLeft = limitProbabilityRight;
        double limitParameterLeft = limitParameterRight;


        // Difference distribution to find point belonging to the distribution to the left
        // of the mode having closest probability to that identified in the part of the
        // distribution to the right of the mode.

        ArrayXd marginalDifferenceLeft = ArrayXd::Zero(max + 1);          
                                                

        // Copy left-hand part of total marginal distributon into separate array.
        // Also copy the corresponding parameter values inro a separate array.
        // Since parameterValues is already sorted in ascending order this can be done quickly.

        marginalDistributionLeft = marginalDistribution.segment(0, max+1);
        parameterValuesLeft = parameterValues.segment(0, max+1);


        // Count number of bins in the distribution to the righti side of the modal value.
        // Define the limit probability on the right-hand distribution 
        // and the corresponding value of the parameter.
        // Then compute the difference of the left-hand distribution with the right limit probability.
        // Finally find the point in the left part having closest probability to the right limit,
        // save its probability and the corresponding value of the parameter.
        // Compute the total probability within the range and repeat the process until 
        // probability >= probability(credibleLevel).

        int stepRight = 0;          
        int min = 0;
        int sampleSize = parameterValues.size();

        while (totalProbability < (credibleLevel/100.) && (max + stepRight < sampleSize))             
        {
            totalProbability = 0.0;


            // Find the probability and the corresponding value of the parameter at the right edge of the interval

            limitProbabilityRight = marginalDistribution(max + stepRight);  
            limitParameterRight = parameterValues(max + stepRight);         

            
            // Find which point in the left part is the closest in probability to that of the right edge.

            marginalDifferenceLeft = (marginalDistributionLeft - limitProbabilityRight).abs();                                      
            limitProbabilityLeft = marginalDifferenceLeft.minCoeff(&min);       


            // Save the point and its probability is the left edge of the search interval

            limitProbabilityLeft = marginalDistributionLeft(min);           
            limitParameterLeft = parameterValuesLeft(min);                 


            // Loop within all the bins in the interval identified in order to compute the total probability

            for (int t = min; t <= max + stepRight; ++t)     
            {
                totalProbability += marginalDistribution(t);
            }


            // Move one bin further to the right

            ++stepRight;
        }
        
        double lowerCredibleLimit = limitParameterLeft;      // Save the left
        double upperCredibleLimit = limitParameterRight;     // and the righ credible limit
           
        parameterEstimates(i,4) = lowerCredibleLimit;
        parameterEstimates(i,5) = upperCredibleLimit;
        
    }   // END for loop over the parameters

    return parameterEstimates;
}











// Results::writeParametersToFile()
//
// PURPOSE: 
//      writes the parameters values from the nested sampling
//      in separate ASCII files having one column format each.
//      
// OUTPUT:
//      void
// 

void Results::writeParametersToFile(string pathPrefix, string outputFileExtension)
{
    ArrayXXd posteriorSample = nestedSampler.getPosteriorSample();
    File::arrayXXdRowsToFiles(posteriorSample, pathPrefix, outputFileExtension);
}












// Results::writeLogLikelihoodToFile()
//
// PURPOSE:
//      writes the log likelihood values from the nested sampling
//      into an ASCII file of one column format. The values are
//      sorted in increasing order.
//
// INPUT:
//      fullPath:   a string variable containing the desired full path to save the output file.
//
// OUTPUT:
//      void
// 

void Results::writeLogLikelihoodToFile(string fullPath)
{
    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
            
    outputFile << "# Posterior sample from nested sampling" << endl;
    outputFile << "# log(Likelihood)" << endl;
    outputFile << scientific << setprecision(9);
    
    ArrayXd logLikelihoodOfPosteriorSample = nestedSampler.getLogLikelihoodOfPosteriorSample();
    File::arrayXdToFile(outputFile, logLikelihoodOfPosteriorSample);
    outputFile.close();
}













// Results::writeLogWeightsToFile()
//
// PURPOSE:
//      writes the log(Weight) = log(dX) values from the nested sampling
//      into an ASCII file of one column format. The values are sorted according to the
//      ascending order in likelihood.
//
// INPUT:
//      fullPath:   a string variable containing the desired full path to save the output file.
//
// OUTPUT:
//      void
// 

void Results::writeLogWeightsToFile(string fullPath)
{
    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
            
    outputFile << "# Posterior sample from nested sampling" << endl;
    outputFile << "# log(Weight) = log(dX)" << endl;
    outputFile << scientific << setprecision(9);
    
    ArrayXd logWeightOfPosteriorSample = nestedSampler.getLogWeightOfPosteriorSample();
    File::arrayXdToFile(outputFile, logWeightOfPosteriorSample);
    outputFile.close();
}













// Results::writeEvidenceInformationToFile()
//
// PURPOSE:
//      writes all the information related to the evidence computed from 
//      the nested sampling into an ASCII file. 
//      Skilling's evidence, evidence error and information gain are also included.
//
// INPUT:
//      fullPath:   a string variable containing the desired full path to save the output file.
//
// OUTPUT:
//      void
//

void Results::writeEvidenceInformationToFile(string fullPath)
{
    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
            
    outputFile << "# Evidence results from nested sampling" << endl;
    outputFile << scientific << setprecision(9);
    outputFile << "# Skilling's log(Evidence)" << setw(40) << "Skilling's Error log(Evidence)" 
    << setw(40) << "Skilling's Information Gain" << endl;
    outputFile << nestedSampler.getLogEvidence() << setw(40) << nestedSampler.getLogEvidenceError() 
    << setw(40) << nestedSampler.getInformationGain() << endl;
    outputFile.close();
} 












// Results::writePosteriorProbabilityToFile()
//
// PURPOSE:
//      writes the posterior probability for the sample obtained from the 
//      nested sampling into an ASCII file of one column format.
//
// INPUT:
//      fullPath:   a string variable containing the desired full path to save the output file.
//
// OUTPUT:
//      void
//
// REMARK:
//      Note that these values are probabilities and not probability densities.
// 

void Results::writePosteriorProbabilityToFile(string fullPath)
{
    ArrayXd posteriorDistribution = posteriorProbability();

    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
            
    outputFile << "# Posterior probability distribution from nested sampling" << endl;
    outputFile << scientific << setprecision(9);
    File::arrayXdToFile(outputFile, posteriorDistribution);
    outputFile.close();

} 









// Results:writeParametersSummaryToToFile()
//
// PURPOSE:
//      Writes the expectation, standard deviation, median and mode values from the 
//      marginalized posterior probability into an ASCII file. 
//      Shortest Bayesian credible intervals (CI) are also included.
//
// INPUT:
//      fullPath:       a string variable containing the desired full path to save the output file.
//      credibleLevel:  a double number providing the desired credible 
//                      level to be computed. Default value corresponds 
//                      to a credible level of 68.27 %.
//      
// OUTPUT:
//
// 

void Results::writeParametersSummaryToFile(string fullPath, const double credibleLevel)
{
    ArrayXXd parameterEstimates = parameterEstimation(credibleLevel);


    // Write output ASCII file

    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
            
    outputFile << "# Summary of Parameter Estimation from nested sampling" << endl;
    outputFile << "# Credible intervals are the shortest credible intervals" << endl; 
    outputFile << "# according to the usual definition" << endl;
    outputFile << "# Credible level: " << fixed << setprecision(2) << credibleLevel << " %" << endl;
    outputFile << "# Column #1: Expectation (I Moment)" << endl;
    outputFile << "# Column #2: Median" << endl;
    outputFile << "# Column #3: Mode" << endl;
    outputFile << "# Column #4: II Moment (Variance if Normal Distribution)" << endl;
    outputFile << "# Column #5: Lower Credible Limit" << endl;
    outputFile << "# Column #6: Upper Credible Limit" << endl;
    outputFile << scientific << setprecision(9);
    File::arrayXXdToFile(outputFile, parameterEstimates);
    outputFile.close();

} 



