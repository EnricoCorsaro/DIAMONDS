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










// Results:parameterEstimation()
//
// PURPOSE:
//      Computes the expectation, median and mode values from the 
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
//      A bidimensional Eigen Array containing all the estimators of the
//      free parameters (one parameter per row).
// 

ArrayXXd Results::parameterEstimation(const double credibleLevel)
{
    int Ndimensions = nestedSampler.posteriorSample.rows();
    int Npoints = nestedSampler.posteriorSample.cols();
    ArrayXd posteriorDistribution = posteriorProbability();
    assert (nestedSampler.posteriorSample.cols() == posteriorDistribution.size());
    ArrayXXd parameterEstimates(Ndimensions, 6);


    // Loop over all free parameters

    for (int i = 0; i < Ndimensions; ++i)
    {
        ArrayXd parameterValues = nestedSampler.posteriorSample.row(i);
        ArrayXd marginalDistribution = posteriorDistribution;


        // Sort elements of array parameterValues in increasing
        // order and sort elements of array marginalDistribution accordingly
        
        Functions::sortElementsDouble(parameterValues, marginalDistribution); 
        

        // Merge existing posterior values that orrespond 
        // to equal parameter values (if any)

        int NduplicateParameterComponents = 0;

        for (int j = 0; j < Npoints - 1; ++j)
        {  
            if (parameterValues(j) != numeric_limits<double>::lowest())
            {
                for (int k = j + 1; k < Npoints; k++)
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


        // Remove bad points and store final values into array copies.
        // Check if bad points are present otherwise skip block.

        if (NduplicateParameterComponents > 0)      
        {
            ArrayXd parameterValuesCopy(Npoints - NduplicateParameterComponents);
            ArrayXd marginalDistributionCopy(Npoints - NduplicateParameterComponents);
        
            int n = 0;

            for (int m = 0; (m < Npoints) && (n < Npoints - NduplicateParameterComponents); m++)
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

        
        // Compute third moment of the sample distribution (including skewness)
        
        double thirdMoment = ((parameterValues - parameterMean).pow(3) * marginalDistribution).sum();
        double skewness = thirdMoment/pow(secondMoment,3.0/2.0);
        

        // Compute optimal bin size for rebinning marginal distribution according to its symmetry properties

        int Nbins = 0;
        int sampleSize = parameterValues.size();
        double binWidth = 0;
        double parameterMaximum = parameterValues.maxCoeff();
        double parameterMinimum = parameterValues.minCoeff();
        
        if (fabs(skewness) < 0.10)
        {
            // The distribution is rather symmetric.
            // Scott's normal reference rule is adopted.
            
            binWidth = 3.5*sqrt(secondMoment)/pow(sampleSize,1.0/3.0);
            Nbins = ceil((parameterMaximum - parameterMinimum)/binWidth);
        }
        else
        {
            // The distribution is far from being symmetric.
            // Doane's formula for non-normal distributions is adopted.

            double constant = sqrt(6.0*(sampleSize - 2.0) / ((sampleSize + 1.0)*(sampleSize + 3.0)));
            double k = 1.0 + log(sampleSize)/log(2) + log(1 + (fabs(skewness)/constant))/log(2);
            Nbins = ceil(k);
            binWidth = (parameterMaximum - parameterMinimum)/k;
        }

        parameterEstimates(i,3) = secondMoment;


        // Rebin marginal distribution for computing credible intervals        

        int beginIndex = 0;
        ArrayXd marginalDistributionRebinned = ArrayXd::Zero(Nbins);
        ArrayXd parameterValuesRebinned(Nbins);


        // Merge marginal distribution for each bin up to the second last bin

        for (int j = 0; j < Nbins-1; ++j)
        {
            parameterValuesRebinned(j) = parameterMinimum + j*binWidth + 0.5*binWidth;
            
            for (int h = beginIndex; h < sampleSize; ++h)
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
            
        for (int h = beginIndex; h < sampleSize; ++h)
        {
            if ((parameterValues(h) >= parameterMinimum + (Nbins-1)*binWidth) && (parameterValues(h) <= parameterMaximum))
            {
                marginalDistributionRebinned(Nbins-1) += marginalDistribution(h);
            }
        }


        // Compute the "shortest" credible intervals (CI)

        marginalDistributionMode = marginalDistributionRebinned.maxCoeff(&max);
        ArrayXd marginalDistributionLeft(max + 1);       // Marginal distribution up to mode value
        ArrayXd parameterValuesLeft(max + 1);            // Parameter range up to mode value
        double limitProbabilityRight = marginalDistributionRebinned(max);
        double limitParameterRight = parameterValuesRebinned(max);
        double limitProbabilityLeft = limitProbabilityRight;
        double limitParameterLeft = limitParameterRight;


        // Difference distribution to find point belonging to the distribution the left of mode 
        // value having closest probability to that identified in the part of the
        // distribution to the right of the mode value.

        ArrayXd marginalDifferenceLeft = marginalDistributionLeft;          
                                                

        // Consider left-hand distribution (up to mode value)
        // Save that part of the parameter values corresponding to the
        // left-hand distribution

        for (int nn = 0; nn <= max; ++nn)
        {
            marginalDistributionLeft(nn) = marginalDistributionRebinned(nn);
            parameterValuesLeft(nn) = parameterValuesRebinned(nn);                                                                                                  }


        // Count number of steps in the distribution to the right of the modal value
        // and stop when probability >= credibleLevel probability.
        // Define the limit probability on the right-hand distribution 
        // and the corresponding value of the parameter.
        // Then compute the difference of the left-hand distribution to the limit probability.
        // Finally find the point having closest probability to the limit and save its 
        // probability and the corresponding value of the parameter.

        int stepRight = 0;          
        int min = 0;

        while (totalProbability < (credibleLevel/100.) && (max + stepRight < Nbins))             
        {

            totalProbability = 0.0;
            limitProbabilityRight = marginalDistributionRebinned(max + stepRight);  
            limitParameterRight = parameterValuesRebinned(max + stepRight);         
            marginalDifferenceLeft = (marginalDistributionLeft - limitProbabilityRight).abs();                                                                                                         
            limitProbabilityLeft = marginalDifferenceLeft.minCoeff(&min);       
            limitProbabilityLeft = marginalDistributionRebinned(min);           
            limitParameterLeft = parameterValuesRebinned(min);                 


            // Loop over the interval identified in order to cumulate the probability and compute the total one

            for (int t = min; t <= max + stepRight; ++t)     
            {
                totalProbability += marginalDistributionRebinned(t);        // Evaluate total probability within the range
            }

            ++stepRight;
        }
        
        double lowerCredibleInterval = parameterMean - limitParameterLeft;      // Compute the left-hand side ...
        double upperCredibleInterval = limitParameterRight - parameterMean;     // and the right-hand side of the credible interval
           
        parameterEstimates(i,4) = lowerCredibleInterval;
        parameterEstimates(i,5) = upperCredibleInterval;

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
    outputFile << "# Column #1: Expectation (I Moment)" << endl;
    outputFile << "# Column #2: Median" << endl;
    outputFile << "# Column #3: Mode" << endl;
    outputFile << "# Column #4: II Moment (Variance if Normal Distribution)" << endl;
    outputFile << "# Column #5: Lower Credible Interval (CI)" << endl;
    outputFile << "# Column #6: Upper Credible Interval (CI)" << endl;
    outputFile << scientific << setprecision(9);
    File::arrayXXdToFile(outputFile, parameterEstimates);
    outputFile.close();

} 



