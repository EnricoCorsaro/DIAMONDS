#include "Results.h"


// Results::Results()
//
// PURPOSE: 
//      Constructor. Sets nested sampler object and output directory path.
//
// INPUT:
//      nestedSampler: a NestedSampler class object used as the container of
//      information to write from.
//

Results::Results(NestedSampler &nestedSampler)
: nestedSampler(nestedSampler)
{
} // END Results::Results()









// Results::~Results()
//
// PURPOSE: 
//      Destructor.
//

Results::~Results()
{

} // END Results::~Results()










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
    int Ndimensions = nestedSampler.posteriorSample.rows();
    assert(Ndimensions > 0);

    // Find out the number of decimal digits that the number of dimensions has
    
    int Ndigits = int(floor(log10(double(Ndimensions)))); 
    
    // Write everything to the output file

    for (int i = 0; i < Ndimensions; i++)
    {
        // Include the dimension serial number with preceding zeros
        
        ostringstream numberString;
        numberString << setfill('0') << setw(Ndigits) << i;
        string fullPath = pathPrefix + numberString.str() + outputFileExtension;
        
        // Open the output file and check for sanity
        
        ofstream outputFile(fullPath.c_str());
        
        if (!outputFile.good())
        {
                cerr << "Error opening output file " << fullPath << endl;
                exit(EXIT_FAILURE);
        }
    
        // Write a header to the output file
        
        outputFile << "# Posterior sample from nested algorithm" << endl;
        outputFile << "# Parameter " + numberString.str() << endl;
        
        // Write all values of this particular parameter in our sample to the output file
        
        outputFile << setiosflags(ios::scientific) << setprecision(9);
        File::oneArrayToFile(outputFile, nestedSampler.posteriorSample.row(i));
        outputFile.close();
    }

} // END Results::writeParametersToFile()









// Results::writeLogLikelihoodToFile()
//
// PURPOSE:
//      writes the log likelihood values from the nested sampling
//      into an ASCII file of one column format. The values are
//      sorted in increasing order.
//
// OUTPUT:
//      void
// 

void Results::writeLogLikelihoodToFile(string fullPath)
{
    ofstream outputFile(fullPath.c_str());
            
    if (!outputFile.good())
    {
        cerr << "Error opening output file " << fullPath << endl;
        exit(EXIT_FAILURE);
    }
            
    outputFile << "# Posterior sample from nested algorithm" << endl;
    outputFile << "# log Likelihood" << endl;
    outputFile << setiosflags(ios::scientific) << setprecision(9);
    File::oneArrayToFile(outputFile, nestedSampler.logLikelihoodOfPosteriorSample);
    outputFile.close();

} // END Results::writeLogLikelihoodToFile()










// Results::writeEvidenceToFile()
//
// PURPOSE:
//      writes the evidence from the nested sampling into an ASCII file. 
//      Evidence error and information Gain are also included.
//
// OUTPUT:
//      void
//
// TODO: rename this function. Current name doesn't cover what it does.
// 

void Results::writeEvidenceToFile(string fullPath)
{
    ofstream outputFile(fullPath.c_str());
            
    if (!outputFile.good())
    {
        cerr << "Error opening output file"  << fullPath << endl;
        exit(EXIT_FAILURE);
    }
            
    outputFile << "# Evidence results from nested algorithm" << endl;
    outputFile << "# log(Evidence)    Error of log(Evidence)    Information Gain" << endl;
    outputFile << setiosflags(ios::scientific) << setprecision(9);
    outputFile << nestedSampler.getLogEvidence() << "    ";
    outputFile << nestedSampler.getLogEvidenceError() << "    ";
    outputFile << nestedSampler.getInformationGain() << endl;
    outputFile.close();

} // END Results::writeEvidenceToFile()









// Results::writePosteriorToFile()
//
// PURPOSE:
//      writes the posterior probability for the sample obtained from the 
//      nested sampling into an ASCII file of one column format. 
//      The values are also stored into the one dimensional Eigen Array
//      posteriorDistribution.
//
// OUTPUT:
//      void
// 

void Results::writePosteriorToFile(string fullPath)
{
    ArrayXd logPosterior = nestedSampler.logWeightOfPosteriorSample - nestedSampler.getLogEvidence();
    posteriorDistribution = logPosterior.exp();

    ofstream outputFile(fullPath.c_str());

    if (!outputFile.good())
    {
        cerr << "Error opening output file"  << fullPath << endl;
        exit(EXIT_FAILURE);
    }
            
    outputFile << "# Posterior probability distribution from nested algorithm" << endl;
    outputFile << scientific << setprecision(9);
    File::oneArrayToFile(outputFile, posteriorDistribution);
    outputFile.close();

} // END Results::writePosteriorToFile()









// Results:writeSummaryStatisticsToFile()
//
// PURPOSE:
//      computes the expectation, median and mode values from the 
//      marginalized posterior probability into an ASCII file. 
//      Shortest Bayesian credible intervals (CI) are also computed.
//      All the values are stored in to the Eigen Array summaryStatistics.
//
// INPUT:
//      credibleLevel: a double number providing the desired credible 
//      level to be computed. Default value correspond to most 
//      used credible level of 68.27 %.
//      
// OUTPUT:
//
// TODO: - Separate the computing of the summary statistics from the writing to the file.
//      
// 

void Results::writeSummaryStatisticsToFile(string fullPath, const double credibleLevel)
{
    int Ndimensions = nestedSampler.posteriorSample.rows();
    int Niterations = nestedSampler.posteriorSample.cols();
    ArrayXd parameterComponent;
    ArrayXd marginalDistribution;

    summaryStatistics.resize(Ndimensions, 5);


    // Loop over all free parameters

    for (int i = 0; i < Ndimensions; i++)
    {
        parameterComponent = nestedSampler.posteriorSample.row(i);
        marginalDistribution = posteriorDistribution;


        // Sort elements of array parameterComponent in increasing
        // order and sort elements of array marginalDistribution accordingly
        
        Functions::sortElements(parameterComponent, marginalDistribution); 
        

        // Marginalize over parameter by merging posterior values 
        // corresponding to equal parameter values

        int NduplicateParameterComponents = 0;

        for (int j = 0; j < Niterations - 1; j++)
        {  
            if (parameterComponent(j) == -DBL_MAX)
                continue;
            else
            {
                for (int k = j + 1; k < Niterations; k++)
                {
                    if (parameterComponent(k) == -DBL_MAX)
                        continue;
                    else
                        if (parameterComponent(j) == parameterComponent(k))
                        {   
                            parameterComponent(k) = -DBL_MAX;        // Set duplicate to bad value (flag)
                            marginalDistribution(j) = marginalDistribution(j) + marginalDistribution(k); // Merge probability values
                            marginalDistribution(k) = 0.0;
                            NduplicateParameterComponents++;
                        }
                }
            }
        }


        // Remove bad points and store final values in array copies

        if (NduplicateParameterComponents > 0) // Check if bad points are present otherwise skip block
        {
            ArrayXd parameterComponentCopy(Niterations - NduplicateParameterComponents);
            ArrayXd marginalDistributionCopy(Niterations - NduplicateParameterComponents);
        
            int n = 0;

            for (int m = 0; (m < Niterations) && (n < (Niterations - NduplicateParameterComponents)); m++)
            {
                if (parameterComponent(m) == -DBL_MAX)
                    continue;
                else
                    if (parameterComponent(m) != -DBL_MAX)
                        {
                            parameterComponentCopy(n) = parameterComponent(m);
                            marginalDistributionCopy(n) = marginalDistribution(m);
                            n++;
                        }
            }


            // Replace original marginal arrays with array copies

            parameterComponent = parameterComponentCopy;
            marginalDistribution = marginalDistributionCopy;
        }

        // Compute the mean value (expectation value)
       
        double meanParameter;
        
        meanParameter = (parameterComponent.cwiseProduct(marginalDistribution)).sum();
        summaryStatistics(i,0) = meanParameter;


        // Compute the median value (value corresponding to 50% of probability)

        double medianParameter;
        double totalProbability = 0.0;
        int k = 0;
        
        while (totalProbability < 0.5)
        {
            medianParameter = parameterComponent(k);
            totalProbability += marginalDistribution(k);
            k++;
        }
        
        summaryStatistics(i,1) = medianParameter;


        // Find the mode value (parameter corresponding to maximum probability value)

        int max = 0;                                    // Subscript corresponding to mode value
        double maximumMarginal;
        double maximumParameter;

        maximumMarginal = marginalDistribution.maxCoeff(&max);
        maximumParameter = parameterComponent(max);
        summaryStatistics(i,2) = maximumParameter;

        
        // Compute the "shortest" credible intervals (CI)

        int stepRight = 1;
        double limitProbabilityRight;
        double limitProbabilityLeft;
        double limitParameterLeft;
        double limitParameterRight;
        ArrayXd marginalDistributionLeft(max + 1);       // Marginal distribution up to mode value
        ArrayXd parameterComponentLeft(max + 1);          // Parameter range up to mode value
        ArrayXd marginalDifferenceLeft;                  // Difference distribution to find point in 
                                                         // in the left-hand distribution having equal
                                                         // probability to that identified in the right-hand part
        for (int nn = 0; nn <= max; nn++)
        {
            marginalDistributionLeft(nn) = marginalDistribution(nn);          // Consider left-hand distribution (up to mode value)
            parameterComponentLeft(nn) = parameterComponent(nn);
        }


        while (totalProbability < (credibleLevel/100.))     // Stop when probability >= credibleLevel 
        {
            totalProbability = 0.0;
            limitProbabilityRight = marginalDistribution(max + stepRight);
            limitParameterRight = parameterComponent(max + stepRight);
            marginalDifferenceLeft = (marginalDistributionLeft - limitProbabilityRight).abs();

            int min = 0;

            limitProbabilityLeft = marginalDifferenceLeft.minCoeff(&min);
            limitProbabilityLeft = marginalDistribution(min);
            limitParameterLeft = parameterComponent(min);

            for (int t = min; t <= (max + stepRight); t++)
            {
                totalProbability += marginalDistribution(t);        // Evaluate total probability within the range
            }

            stepRight++;
        }

        double lowerCredibleInterval;
        double upperCredibleInterval;

        lowerCredibleInterval = maximumParameter - limitParameterLeft;
        upperCredibleInterval = limitParameterRight - maximumParameter;
           
        summaryStatistics(i,3) = lowerCredibleInterval;
        summaryStatistics(i,4) = upperCredibleInterval;

    }


    // Write output ASCII file

    ofstream outputFile(fullPath.c_str());

    if (!outputFile.good())
    {
        cerr << "Error opening output file" << endl;
        exit(EXIT_FAILURE);
    }
            
    outputFile << "# Summary statistics from MultiNest algorithm" << endl;
    outputFile << "# Credible intervals are the shortest credible intervals" << endl; 
    outputFile << "# according to the usual definition" << endl;
    outputFile << "# Credible level: " << fixed << setprecision(2) << credibleLevel << " %" << endl;
    outputFile << "# Column #1: Expectation" << endl;
    outputFile << "# Column #2: Median" << endl;
    outputFile << "# Column #3: Mode" << endl;
    outputFile << "# Column #4: Lower Credible Interval (CI)" << endl;
    outputFile << "# Column #5: Upper Credible Interval (CI)" << endl;
    outputFile << fixed << setprecision(12);
    File::arrayToFile(outputFile, summaryStatistics);
    outputFile.close();

} // END Results::writeSummaryStatisticsToFile()


















// Results::getPosteriorDistribution()
//
// PURPOSE:
//      Gets private data member posteriorDistribution.
//
// OUTPUT:
//      An Eigen Array containing the values of the posterior
//      distribution.
//
// TODO: fixme: when this member function is called before any other member function,
//              the result is undefined
// 

ArrayXd Results::getPosteriorDistribution()
{
    return posteriorDistribution;

} // END Results::getPosteriorDistribution()











// Results::getSummaryStatistics()
//
// PURPOSE:
//      Gets private data member summaryStatistics.
//
// OUTPUT:
//      An Eigen Array containing all the values obtained
//      from the inference analysis of the posterior probability
//      distribution.
//
// TODO: - fixme: when this member function is called before any other member function,
//                the result is undefined
//       - fixme: returning an array for which the user has to decipher what quantity is 
//                located in which array element is not the best practice. Find a better
//                solution
//       
// 

ArrayXXd Results::getSummaryStatistics()
{
    return summaryStatistics;

} // END Results::getSummaryStatistics()


