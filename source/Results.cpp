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
//      parameterNumber:            an integer containing the number of the free parameter to save, together with
//                                  its marginal distribution.
//      
// OUTPUT:
//      void.
// 

void Results::writeMarginalDistributionToFile(const int parameterNumber)
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
    string fileName = "marginalDistribution";
    string fullPath = nestedSampler.getOutputPathPrefix() + fileName + numberString.str() + ".txt";
    
    
    // Write the two-columns array in an ASCII file

    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
    outputFile << "# Marginal distribution of cubic-spline interpolated points from nested sampling." << endl;
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
//      skewness:                   the skewness of the sample distribution
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

ArrayXd Results::computeCredibleLimits(const double credibleLevel, const double skewness, const int NinterpolationsPerBin)
{
    // Interpolate rebinned marginal distribution by using a NinterpolationsPerBin times finer grid. 
    // This allows for a better computation of the credible intervals. Note that while the rebinned 
    // marginal distribution has still sum equal to 1, the interpolated marginal distribution must 
    // be renormalized after interpolation. This will scale down the amplitudes of a factor of the order of NinterpolationsPerBin. 
    // Such a rescaling has no effect on the computation of the credible intervals,
    // which instead depend upon the relative variation from point to point in the distribution.

    int Nbins = parameterValuesRebinned.size();
    int Ninterpolations = Nbins*NinterpolationsPerBin;


    // Take the average bin width of the rebinned data for computing the bin width of the interpolated data. This information is more
    // realistic as it better resembles the original sampling from the nesting process.

    double binWidth = parameterValuesRebinned(1) - parameterValuesRebinned(0);
    double interpolatedBinWidth = binWidth/NinterpolationsPerBin;
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
    int NbinsLeft = max + 1;
    int NbinsRight = Ninterpolations - NbinsLeft;
    ArrayXd marginalDistributionLeft(NbinsLeft);        // Marginal distribution up to modal value (included)
    ArrayXd parameterValuesLeft(NbinsLeft);             // Parameter range up to modal value (included)
    ArrayXd marginalDistributionRight(NbinsRight);      // Marginal distribution beyond modal value 
    ArrayXd parameterValuesRight(NbinsRight);           // Parameter range beyond modal value
    double limitProbabilityRight = marginalDistributionInterpolated(max);
    double limitParameterRight = parameterValuesInterpolated(max);
    double limitProbabilityLeft = limitProbabilityRight;
    double limitParameterLeft = limitParameterRight;


    // Difference distribution to find point belonging to the distribution to the left (right)
    // of the mode having closest probability to that identified in the part of the
    // distribution to the right (left) of the mode.

    ArrayXd marginalDifferenceLeft = ArrayXd::Zero(NbinsLeft);          
    ArrayXd marginalDifferenceRight = ArrayXd::Zero(NbinsRight);          
                                                

    // Copy left- and right-hand part of total marginal distributon into separate arrays.
    // Also copy the corresponding parameter values inro separate arrays.
    // Since parameterValues is already sorted in ascending order this can be done quickly.

    marginalDistributionLeft = marginalDistributionInterpolated.segment(0, NbinsLeft);
    marginalDistributionRight = marginalDistributionInterpolated.segment(NbinsLeft, NbinsRight);
    parameterValuesLeft = parameterValuesInterpolated.segment(0, NbinsLeft);
    parameterValuesRight = parameterValuesInterpolated.segment(NbinsLeft, NbinsRight);


    // Count number of steps (bins) in the distribution to the right (left) side of its modal value.
    // Define the limit probability on the right-(left-)hand distribution and the corresponding 
    // value of the parameter. Then compute the difference of the left-(right-)hand distribution 
    // with the right (left) limit probability. Finally find the point in the left (right) part 
    // having closest probability to the right (left) limit, save its probability and the corresponding 
    // value of the parameter. Cumulate the probability within the range and repeat the process until 
    // probability >= probability(credibleLevel). The skenwess is adopted to select the direction 
    // of the iterative process for a more optimal computation of the shortest CI. 
    // For skewness < 0, the distribution is asymetric to the left, hence the step is towards the left side.
    // For skewness > 0, the distribution is asymetric to the right, hence the step is towards the right side.

    int min = 0;
    double credibleLevelFraction = credibleLevel/100.;
    double totalProbability = 0.0;

    if (skewness >= 0.0)
    {
        int stepRight = 0;          
        
        while (totalProbability < (credibleLevelFraction) && (NbinsLeft + stepRight < Ninterpolations))             
        {
            totalProbability = 0.0;


            // Find the probability and the corresponding value of the parameter at the right edge of the interval

            limitProbabilityRight = marginalDistributionInterpolated(NbinsLeft + stepRight);  
            limitParameterRight = parameterValuesInterpolated(NbinsLeft + stepRight);         

            
            // Find which point in the left part is the closest in probability to that of the right edge.

            marginalDifferenceLeft = (marginalDistributionLeft - limitProbabilityRight).abs();                                      
            limitProbabilityLeft = marginalDifferenceLeft.minCoeff(&min);       


            // Save the point and its probability as the left edge of the searching interval

            limitProbabilityLeft = marginalDistributionLeft(min);
            limitParameterLeft = parameterValuesLeft(min);


            // Cumulate the probability within the identified interval

            int intervalSize = NbinsLeft + stepRight - min;
            totalProbability = marginalDistributionInterpolated.segment(min, intervalSize).sum();


            // Move one step further to the right

            ++stepRight;
        }
    }
    else
    {
        int stepLeft = 0;          
        
        while (totalProbability < (credibleLevelFraction) && (stepLeft < NbinsLeft))             
        {
            totalProbability = 0.0;


            // Find the probability and the corresponding value of the parameter at the left edge of the interval

            limitProbabilityLeft = marginalDistributionInterpolated(NbinsLeft - stepLeft - 1);  
            limitParameterLeft = parameterValuesInterpolated(NbinsLeft - stepLeft - 1);         

            
            // Find which point in the right part is the closest in probability to that of the left edge.

            marginalDifferenceRight = (marginalDistributionRight - limitProbabilityLeft).abs();                                      
            limitProbabilityRight = marginalDifferenceRight.minCoeff(&min);       


            // Save the point and its probability as the right edge of the searching interval

            limitProbabilityRight = marginalDistributionRight(min);
            limitParameterRight = parameterValuesRight(min);


            // Cumulate the probability within the identified interval

            int intervalSize = min + stepLeft;
            totalProbability = marginalDistributionInterpolated.segment(NbinsLeft - stepLeft, intervalSize).sum();


            // Move one step further to the left

            ++stepLeft;
        }
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
//      credibleLevel:              a double number providing the desired credible 
//                                  level to be computed. Default value corresponds to
//                                  credible level of 68.27 %.
//      writeMarginalDistribution:  a boolean variable specifying whether the marginal distribution for
//                                  each parameter has to be written in an output file.
//      
// OUTPUT:
//      A bidimensional Eigen Array containing all the estimators of the
//      free parameters (one parameter per row). The output column order is:
//      (1) Mean (I Moment)
//      (2) Median
//      (3) Mode
//      (4) II Moment
//      (5) Lower CL
//      (6) Upper CL
// 

ArrayXXd Results::parameterEstimation(double credibleLevel, bool writeMarginalDistribution)
{
    int Ndimensions = nestedSampler.getPosteriorSample().rows();
    ArrayXd posteriorDistribution = posteriorProbability();
    
    int sampleSize = posteriorDistribution.size();
    assert(nestedSampler.getPosteriorSample().cols() == sampleSize);
    ArrayXXd parameterEstimates(Ndimensions, 7);

    parameterValues.resize(sampleSize);
    marginalDistribution.resize(sampleSize);


    // Loop over all free parameters

    for (int i = 0; i < Ndimensions; ++i)
    {
        // Take the information corresponding to the current parameter

        parameterValues = nestedSampler.getPosteriorSample().row(i);
        marginalDistribution = posteriorDistribution;


        // Sort elements of array parameterValues in ascending
        // order and sort elements of array marginalDistribution accordingly.
        // Use a mergesort algorithm for more efficiency with a large number of array elements (O(n*log n)).
        
        Functions::topDownMergeSort(parameterValues, marginalDistribution); 
        

        // Compute the mean value (expectation value, or first moment) of the sample distribution
       
        double parameterMean = (parameterValues * marginalDistribution).sum();
        parameterEstimates(i,0) = parameterMean;


        // Compute second moment of the sample distribution (i.e. the variance in case of a normal distribution)

        double secondMoment = ((parameterValues - parameterMean).pow(2) * marginalDistribution).sum();


        // Compute the median value (value containing the 50% of total probability)

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

       
        // Compute optimal bin size for rebinning marginal distribution according to its symmetry properties
        // Scott's normal reference rule is adopted (most efficient for Gaussian-shaped distributions)

        double binWidth = 0;
        double parameterMaximum = parameterValues.maxCoeff();
        double parameterMinimum = parameterValues.minCoeff();
        
        binWidth = 3.5*sqrt(secondMoment)/pow(sampleSize,1.0/3.0);
        int Nbins = floor((parameterMaximum - parameterMinimum)/binWidth) - 1;

        if (Nbins > 1000)
        {
            // If it happens that the number of bins is too large, reduce it to a fixed value

            Nbins = 1000;
            binWidth = (parameterMaximum - parameterMinimum)/(Nbins*1.0);
        }


        // Rebin marginal distribution. This allows to get rid of the
        // multiple spikes appearing in the final nested sampling posterior. One can then operate on a smoother
        // shape which allows to derive more reliable credible intervals.
        // For this purpose simply cumulate initial marginal distribution values within each bin.
        // This can be done this way because marginal distribution values are probabilities.

        double parameterStart = 0.0;
        double parameterEnd = 0.0;
        parameterValuesRebinned.resize(Nbins);
        parameterValuesRebinned.setZero();
        marginalDistributionRebinned.resize(Nbins);
        marginalDistributionRebinned.setZero();

        int Nshifts = 20;                           // Total number of shifts for the starting point of the rebinning
        double shiftWidth = binWidth/Nshifts;       // Width of the shift bin
        ArrayXd parameterValuesRebinnedPerShift(Nbins);
        parameterValuesRebinnedPerShift.setZero();
        ArrayXd marginalDistributionRebinnedPerShift(Nbins);
        marginalDistributionRebinnedPerShift.setZero();


        // Do the merging
        // First loop over the offset (shift) for the starting point of the rebinning.
        // The larger the number of shifts, the better the averaged result.

        for (int k = 0; k < Nshifts; ++k)
        {
            int cumulatedBinSize = 0;
            
            // Now loop over the different bins for collecting marginal probability in each of them.

            for (int j = 0; j < Nbins; ++j)
            {
                // Set the left edge of the selected bin

                parameterStart = parameterMinimum + j*binWidth + k*shiftWidth;

                // Ensure the right edge is not exceeding the right array boundary
                
                if (j < (Nbins - 1)) 
                    parameterEnd = parameterMinimum + (j+1)*binWidth + k*shiftWidth;
                else
                    parameterEnd = parameterMaximum;

                
                // Find the number of array elements belonging to the selected bin and take as parameter value the mid point

                int binSize = Functions::countArrayIndicesWithinBoundaries(parameterValues, parameterStart, parameterEnd);
                parameterValuesRebinnedPerShift(j) = (parameterStart + parameterEnd)/2.0;

                if (binSize > 0)
                {
                    // At least one point is found in this bin, hence cumulate the marginal distribution values 
                    // falling inside the selected bin
                    
                    marginalDistributionRebinnedPerShift(j) = marginalDistribution.segment(cumulatedBinSize, binSize).sum();
                    cumulatedBinSize += binSize;
                }
                else
                {
                    // No points are found in this bin, hence set marginal probability to zero

                    marginalDistributionRebinnedPerShift(j) = 0.0;
                }
                
            }
        
            
            // Cumulate the values of the rebinning into a total array

            parameterValuesRebinned += parameterValuesRebinnedPerShift;
            marginalDistributionRebinned += marginalDistributionRebinnedPerShift;
        }


        // Average all the rebinnings done by the total number of shifts adopted
        
        parameterValuesRebinned /= Nshifts;
        marginalDistributionRebinned /= Nshifts;


        // Find the mode value (parameter corresponding to maximum probability value of the rebinned distribution)

        int max = 0;                                    // Subscript corresponding to mode value
        marginalDistributionMode = marginalDistributionRebinned.maxCoeff(&max);
        double parameterMode = parameterValuesRebinned(max);
        parameterEstimates(i,2) = parameterMode;

        
        // Save the second moment of the distribution
       
        parameterEstimates(i,3) = secondMoment;


        // Compute third moment and skewness of the sample distribution
        
        double thirdMoment = ((parameterValues - parameterMean).pow(3) * marginalDistribution).sum();
        double skewness = thirdMoment/pow(secondMoment,3.0/2.0);

        
        // Compute shortest credible intervals (CI) and save their corresponding limiting values (credible limits)

        ArrayXd credibleLimits(2);
        credibleLimits = computeCredibleLimits(credibleLevel, skewness);
        
        parameterEstimates(i,4) = credibleLimits(0);
        parameterEstimates(i,5) = credibleLimits(1);


        // Save the skewness of the distribution

        parameterEstimates(i,6) = skewness;

        
        // If required, save the interpolated marginal distribution in an output file

        if (writeMarginalDistribution)
        {
            writeMarginalDistributionToFile(i);
        }

        
    }   // END for loop over the parameters

    return parameterEstimates;
}












// Results::writeParametersToFile()
//
// PURPOSE: 
//      writes the parameters values from the nested sampling
//      in separate ASCII files having one column format each.
//      
// INPUT:
//      fileName:   a string variable containing the file name of the output file to be saved.
//
// OUTPUT:
//      void
// 

void Results::writeParametersToFile(string fileName, string outputFileExtension)
{
    string pathPrefix = nestedSampler.getOutputPathPrefix() + fileName;
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
//      fileName:   a string variable containing the file name of the output file to be saved.
//
// OUTPUT:
//      void
// 

void Results::writeLogLikelihoodToFile(string fileName)
{
    string fullPath = nestedSampler.getOutputPathPrefix() + fileName;

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
//      fileName:   a string variable containing the file name of the output file to be saved.
//
// OUTPUT:
//      void
// 

void Results::writeLogWeightsToFile(string fileName)
{
    string fullPath = nestedSampler.getOutputPathPrefix() + fileName;

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
//      fileName:   a string variable containing the file name of the output file to be saved.
//
// OUTPUT:
//      void
//

void Results::writeEvidenceInformationToFile(string fileName)
{
    string fullPath = nestedSampler.getOutputPathPrefix() + fileName;

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
//      fileName:   a string variable containing the file name of the output file to be saved.
//
// OUTPUT:
//      void
//
// REMARK:
//      Note that these values are probabilities and not probability densities.
// 

void Results::writePosteriorProbabilityToFile(string fileName)
{
    ArrayXd posteriorDistribution = posteriorProbability();
    string fullPath = nestedSampler.getOutputPathPrefix() + fileName;
    
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
//      fileName:                   a string variable containing the filename of the output file.
//      credibleLevel:              a double number providing the desired credible 
//                                  level to be computed. Default value corresponds 
//                                  to a credible level of 68.27 %.
//      writeMarginalDistribution:  a boolean variable specifying whether the marginal distribution for
//                                  each parameter has to be written in an output file.
//      
// OUTPUT:
//      void.
// 

void Results::writeParametersSummaryToFile(string fileName, const double credibleLevel, const bool writeMarginalDistribution)
{
    // Compute estimators for all the free parameters

    ArrayXXd parameterEstimates = parameterEstimation(credibleLevel, writeMarginalDistribution);


    // Write output ASCII file

    string fullPath = nestedSampler.getOutputPathPrefix() + fileName;
    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
            
    outputFile << "# Summary of Parameter Estimation from nested sampling" << endl;
    outputFile << "# Credible intervals are the shortest credible intervals" << endl; 
    outputFile << "# according to the usual definition" << endl;
    outputFile << "# Credible level: " << fixed << setprecision(2) << credibleLevel << " %" << endl;
    outputFile << "# Column #1: I Moment (Mean)" << endl;
    outputFile << "# Column #2: Median" << endl;
    outputFile << "# Column #3: Mode" << endl;
    outputFile << "# Column #4: II Moment (Variance if Normal distribution)" << endl;
    outputFile << "# Column #5: Lower Credible Limit" << endl;
    outputFile << "# Column #6: Upper Credible Limit" << endl;
    outputFile << "# Column #7: Skewness (Asymmetry of the distribution, -1 to the left, +1 to the right, 0 if symmetric)" << endl;
    outputFile << scientific << setprecision(9);
    File::arrayXXdToFile(outputFile, parameterEstimates);
    outputFile.close();

} 




