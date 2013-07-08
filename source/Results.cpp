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
//      Saves the posterior probability for the sample obtained from the 
//      nested sampling into a one dimensional Eigen Array. 
//
// OUTPUT:
//      An Eigen Array containing the values of the posterior probability
//      sorted according to the nesting algorithm process.
// 
// REMARK:
//      Values are real probabilities (and not probability densities).
//

ArrayXd Results::posteriorProbability()
{
    ArrayXd logPosteriorDistribution = nestedSampler.logWeightOfPosteriorSample - nestedSampler.getLogEvidence();
    ArrayXd posteriorDistribution = logPosteriorDistribution.exp();
    
    return posteriorDistribution;
}










// Results:writeParameterEstimationToFile()
//
// PURPOSE:
//      computes the expectation, median and mode values from the 
//      marginalized posterior probability. 
//      Shortest Bayesian credible intervals (CI) are also computed.
//      All the values are stored in a bidimensional Eigen Array.
//
// INPUT:
//      credibleLevel:  a double number providing the desired credible 
//                      level to be computed. Default value corresponds to
//                      credible level of 68.27 %.
//      
// OUTPUT:
//      A bidimensional Eigen Array containing all the estimates of the
//      free parameters (one parameter per row).
// 

ArrayXXd Results::parameterEstimation(const double credibleLevel)
{
    int Ndimensions = nestedSampler.posteriorSample.rows();
    int Niterations = nestedSampler.posteriorSample.cols();
    ArrayXd posteriorDistribution = posteriorProbability();
    ArrayXXd parameterEstimates(Ndimensions, 6);


    // Loop over all free parameters

    for (int i = 0; i < Ndimensions; ++i)
    {
        ArrayXd parameterValues = nestedSampler.posteriorSample.row(i);
        ArrayXd marginalDistribution = posteriorDistribution;


        // Sort elements of array parameterValues in increasing
        // order and sort elements of array marginalDistribution accordingly
        
        Functions::sortElementsDouble(parameterValues, marginalDistribution); 
        

        // Merge existing posterior values 
        // corresponding to equal parameter values (if any)

        int NduplicateParameterComponents = 0;

        for (int j = 0; j < Niterations - 1; ++j)
        {  
            if (parameterValues(j) != numeric_limits<double>::lowest())
            {
                for (int k = j + 1; k < Niterations; k++)
                {
                    if (parameterValues(k) == numeric_limits<double>::lowest())
                        continue;
                    else
                        if (parameterValues(j) == parameterValues(k))
                        {   
                            parameterValues(k) = numeric_limits<double>::lowest();        // Set duplicate to bad value (flag)
                            marginalDistribution(j) = marginalDistribution(j) + marginalDistribution(k); // Merge probability values
                            marginalDistribution(k) = 0.0;
                            NduplicateParameterComponents++;
                        }
                }
            }
            else
                continue;
        }


        // Remove bad points and store final values into array copies

        if (NduplicateParameterComponents > 0)      // Check if bad points are present otherwise skip block
        {
            ArrayXd parameterValuesCopy(Niterations - NduplicateParameterComponents);
            ArrayXd marginalDistributionCopy(Niterations - NduplicateParameterComponents);
        
            int n = 0;

            for (int m = 0; (m < Niterations) && (n < (Niterations - NduplicateParameterComponents)); m++)
            {
                if (parameterValues(m) == numeric_limits<double>::lowest())
                    continue;
                else
                    if (parameterValues(m) != numeric_limits<double>::lowest())
                        {
                            parameterValuesCopy(n) = parameterValues(m);
                            marginalDistributionCopy(n) = marginalDistribution(m);
                            n++;
                        }
            }


            // Replace original marginal arrays with array copies

            parameterValues = parameterValuesCopy;
            marginalDistribution = marginalDistributionCopy;
        }


        // Compute the mean value (expectation value)
       
        double parameterMean = (parameterValues * marginalDistribution).sum();
        parameterEstimates(i,0) = parameterMean;


        // Compute standard deviation of the sample distribution

        double parameterStandardDeviation = sqrt((pow(parameterValues - parameterMean, 2) * marginalDistribution).sum());

        
        // Compute the median value (value corresponding to 50% of probability)

        double totalProbability = 0.0;
        int k = 0;
        
        while (totalProbability < 0.50)
        {
            double parameterMedian = parameterValues(k);
            totalProbability += marginalDistribution(k);
            k++;
        }
        
        parameterEstimates(i,1) = parameterMedian;


        // Find the mode value (parameter corresponding to maximum probability value)

        int max = 0;                                    // Subscript corresponding to mode value
        double marginalDistributionMode = marginalDistribution.maxCoeff(&max);
        double parameterMode = parameterValues(max);
        parameterEstimates(i,2) = parameterMode;

        
        // Compute optimal bin size for rebinning marginal distribution according to Scott's normal reference rule

        double binWidth = 3.5*parameterStandardDeviation/pow(parameterValues.size(),1.0/3.0);
        parameterEstimates(i,3) = parameterStandardDeviation;


        // Rebin marginal distribution for computing credible intervals        

        double parameterMaximum = parameterValues.maxCoeff();
        double parameterMinimum = parameterValues.minCoeff();
        int Nbins = ceil((parameterMaximum - parameterMinimum)/binWidth);
        int beginIndex = 0;
        ArrayXd marginalDistributionRebinned = ArrayXd::Zero(Nbins);
        ArrayXd parameterValuesRebinned(Nbins);


        // Merge marginal distribution for each bin up to the second last bin

        for (int j = 0; j < Nbins-1; ++j)
        {
            parameterValuesRebinned(j) = parameterMinimum + j*binWidth + 0.5*binWidth;
            
            for (int h = beginIndex; h < parameterValues.size(); ++h)
            {
                if ((parameterValues(h) >= parameterMinimum + j*binWidth) && (parameterValues(h) < parameterMinimum + (j+1)*binWidth))
                {
                    marginalDistributionRebinned(j) += marginalDistribution(h);
                    beginIndex = h;         // Save last useful index for the selected bin
                }
            }
        }


        // Merge marginal distribution in the last bin

        parameterValuesRebinned(Nbins-1) = parameterMinimum + (Nbins-1)*binWidth + 0.5*(parameterMaximum - parameterMinimum - (Nbins-1)*binWidth);
            
        for (int h = beginIndex; h < parameterValues.size(); ++h)
        {
            if ((parameterValues(h) >= parameterMinimum + (Nbins-1)*binWidth) && (parameterValues(h) <= parameterMaximum))
            {
                marginalDistributionRebinned(Nbins-1) += marginalDistribution(h);
            }
        }


        // Compute the "shortest" credible intervals (CI)

        int stepRight = 1;          // Count number of steps in the distribution to the right of the modal value
        marginalDistributionMode = marginalDistributionRebinned.maxCoeff(&max);
        ArrayXd marginalDistributionLeft(max + 1);       // Marginal distribution up to mode value
        ArrayXd parameterValuesLeft(max + 1);            // Parameter range up to mode value
        
        for (int nn = 0; nn <= max; ++nn)
        {
            marginalDistributionLeft(nn) = marginalDistributionRebinned(nn);          // Consider left-hand distribution (up to mode value)
            parameterValuesLeft(nn) = parameterValuesRebinned(nn);
        }

        while (totalProbability < (credibleLevel/100.))     // Stop when probability >= credibleLevel probability
        {
            totalProbability = 0.0;
            double limitProbabilityRight = marginalDistributionRebinned(max + stepRight);
            double limitParameterRight = parameterValuesRebinned(max + stepRight);
            ArrayXd marginalDifferenceLeft = (marginalDistributionLeft - limitProbabilityRight).abs();        // Difference distribution to find point 
                                                                                                              // in the distribution to the left of mode 
                                                                                                              // value having closest probability to 
                                                                                                              // that identified in the part of the
                                                                                                              // distribution to the right of the mode 
                                                                                                              // value.

            int min = 0;
            double limitProbabilityLeft = marginalDifferenceLeft.minCoeff(&min);
            limitProbabilityLeft = marginalDistributionRebinned(min);
            double limitParameterLeft = parameterValuesRebinned(min);

            for (int t = min; (t <= (max + stepRight)) && (t < Nbins); ++t)
            {
                totalProbability += marginalDistributionRebinned(t);        // Evaluate total probability within the range
            }

            ++stepRight;
        }
        
        double lowerCredibleInterval = parameterMean - limitParameterLeft;
        double upperCredibleInterval = limitParameterRight - parameterMean;
           
        parameterEstimates(i,3) = lowerCredibleInterval;
        parameterEstimates(i,4) = upperCredibleInterval;

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
    File::arrayXXdRowsToFiles(nestedSampler.posteriorSample, pathPrefix, outputFileExtension);
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
    File::arrayXdToFile(outputFile, nestedSampler.logLikelihoodOfPosteriorSample);
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
//                      level to be computed. Default value correspond to most 
//                      used credible level of 68.27 %.
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
    outputFile << "# Column #1: Expectation" << endl;
    outputFile << "# Column #2: Median" << endl;
    outputFile << "# Column #3: Mode" << endl;
    outputFile << "# Column #4: Standard Deviation" << endl;
    outputFile << "# Column #5: Lower Credible Interval (CI)" << endl;
    outputFile << "# Column #6: Upper Credible Interval (CI)" << endl;
    outputFile << scientific << setprecision(9);
    File::arrayXXdToFile(outputFile, parameterEstimates);
    outputFile.close();

} 



