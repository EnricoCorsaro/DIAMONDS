#include "Results.h"


// Results::Results()
//
// PURPOSE: 
//      Constructor. Sets nested sampler object and output directory path.
//
// INPUT:
//      nestedSampler: a NestedSampler class object used as the container of
//      information to print from.
//      outputDirectory: a pointer to an array of chars containing the desidered
//      path name where output files are to be printed
//

Results::Results(NestedSampler &nestedSampler, const char *inputDirectory, const char *inputFilename, const char *outputDirectory)
: outputDirectory(outputDirectory),
  nestedSampler(nestedSampler)
{
    int length;

    dataFilename = inputFilename;
    length = dataFilename.length();
    dataFilename = dataFilename.substr(0, length-4) + "-";

} // END Results::Results()









// Results::~Results()
//
// PURPOSE: 
//      Destructor.
//

Results::~Results()
{

} // END Results::~Results()










// Results::printParameters()
//
// PURPOSE: 
//      Prints the parameters values from the nested sampling
//      in separate ASCII files having one column format each.
//      
// OUTPUT:
//      void
// 

void Results::printParameters()
{
    string prefix1 = "parameter00";
    string prefix2 = "parameter0";
    string prefix3 = "parameter";
    string postfix = ".txt";
    string filenameString;
    const char *filename = 0;
    int Ndimensions;
    ofstream outputFile;

    Ndimensions = nestedSampler.posteriorSample.rows();
    assert(Ndimensions > 0);

    if (Ndimensions < 10)
    {
        for (int i = 0; i < Ndimensions; i++)
        {
            filenameString = outputDirectory + dataFilename + prefix1 + to_string(i) + postfix;
            filename = filenameString.data();
            outputFile.open(filename);
            
            if (!outputFile.good())
            {
                    cerr << "Error opening output file" << endl;
                    exit(EXIT_FAILURE);
            }
            
            outputFile << "# Posterior sample from nested algorithm" << endl;
            outputFile << "# Parameter 00"+to_string(i) << endl;
            outputFile << setiosflags(ios::fixed) << setprecision(12);
            File::oneArrayToFile(outputFile, nestedSampler.posteriorSample.row(i));
            outputFile.close();
        }
    }
    else
        if ((Ndimensions >= 10) && (Ndimensions < 100))
        {
            for (int i = 0; i < Ndimensions; i++)
            {
                filenameString = outputDirectory + dataFilename + prefix2 + to_string(i) + postfix;
                filename = filenameString.data();
                outputFile.open(filename);

                if (!outputFile.good())
                {
                    cerr << "Error opening output file" << endl;
                    exit(EXIT_FAILURE);
                }

                outputFile << "# Posterior sample from nested algorithm" << endl;
                outputFile << "# Parameter 0"+to_string(i) << endl;
                outputFile << setiosflags(ios::fixed) << setprecision(12);
                File::oneArrayToFile(outputFile, nestedSampler.posteriorSample.row(i));
                outputFile.close();
            } 
        }
    else
        if (Ndimensions >= 100)
        {
            for (int i = 0; i < Ndimensions; i++)
            {
                filenameString = outputDirectory + dataFilename + prefix3 + to_string(i) + postfix;
                filename = filenameString.data();
                outputFile.open(filename);
                
                if (!outputFile.good())
                {
                    cerr << "Error opening output file" << endl;
                    exit(EXIT_FAILURE);
                }
                
                outputFile << "# Posterior sample from nested algorithm" << endl;
                outputFile << "# Parameter "+to_string(i) << endl;
                outputFile << setiosflags(ios::fixed) << setprecision(12);
                File::oneArrayToFile(outputFile, nestedSampler.posteriorSample.row(i));
                outputFile.close();
            }
        }
} // END Results::printParameters()









// Results::printLogLikelihood()
//
// PURPOSE:
//      Prints the log likelihood values from the nested sampling
//      into an ASCII file of one column format. The values are
//      sorted in increasing order.
//
// OUTPUT:
//      void
// 

void Results::printLogLikelihood()
{
    string name = "loglikelihood.txt";
    string filenameString;
    const char *filename = 0;

    filenameString = outputDirectory + dataFilename + name;
    filename = filenameString.data();
    ofstream outputFile(filename);
            
    if (!outputFile.good())
    {
        cerr << "Error opening output file" << endl;
        exit(EXIT_FAILURE);
    }
            
    outputFile << "# Posterior sample from nested algorithm" << endl;
    outputFile << "# log Likelihood" << endl;
    outputFile << setiosflags(ios::fixed) << setprecision(12);
    File::oneArrayToFile(outputFile, nestedSampler.logLikelihoodOfPosteriorSample);
    outputFile.close();

} // END Results::printLogLikelihood()










// Results::printEvidence()
//
// PURPOSE:
//      Prints the evidence from the nested sampling
//      into an ASCII file. Evidence error
//      and information Gain are also included.
//
// OUTPUT:
//      void
// 

void Results::printEvidence()
{
    string name = "evidence.txt";
    string filenameString;
    const char *filename = 0;

    filenameString = outputDirectory + dataFilename + name;
    filename = filenameString.data();
    ofstream outputFile(filename);
            
    if (!outputFile.good())
    {
        cerr << "Error opening output file" << endl;
        exit(EXIT_FAILURE);
    }
            
    outputFile << "# Evidence results from nested algorithm" << endl;
    outputFile << "# log Evidence" << endl;
    outputFile << setiosflags(ios::fixed) << setprecision(12);
    outputFile << nestedSampler.getLogEvidence() << endl;
    outputFile << "# ------------------------ #" << endl;
    outputFile << "# log Evidence Error" << endl;
    outputFile << nestedSampler.getLogEvidenceError() << endl;
    outputFile << "# ------------------------ #" << endl;
    outputFile << "# Information Gain" << endl;
    outputFile << nestedSampler.getInformationGain() << endl;
    outputFile.close();

} // END Results::printEvidence()









// Results::printPosterior()
//
// PURPOSE:
//      Prints the posterior probability for the sample 
//      obtained from the nested sampling into an ASCII file of one column format. 
//      The values are also stored into the one dimensional Eigen Array
//      posteriorDistribution.
//
// OUTPUT:
//      void
// 

void Results::printPosterior()
{
    ArrayXd logPosterior;
    string name = "posterior.txt";
    string filenameString;
    const char *filename = 0;

    logPosterior = nestedSampler.logWeightOfPosteriorSample 
    - nestedSampler.getLogEvidence();
    posteriorDistribution = logPosterior.exp();

    filenameString = outputDirectory + dataFilename + name;
    filename = filenameString.data();
    ofstream outputFile(filename);

    if (!outputFile.good())
    {
        cerr << "Error opening output file" << endl;
        exit(EXIT_FAILURE);
    }
            
    outputFile << "# Posterior probability distribution from nested algorithm" << endl;
    outputFile << scientific << setprecision(12);
    File::oneArrayToFile(outputFile, posteriorDistribution);
    outputFile.close();

} // END Results::printPosterior()










// Results::printInference()
//
// PURPOSE:
//      Prints the expectation, median and mode values from the 
//      marginalized posterior probability into an ASCII file. 
//      Shortest Bayesian credible intervals (CI) are also computed.
//      All the values are stored in to the Eigen Array inferenceResults.

// INPUT:
//      credibleLevel: a double number providing the desired credible 
//      level to be computed. Default value correspond to most 
//      used credible level of 68.27 %.
//      
// OUTPUT:
//      void
// 

void Results::printInference(const double credibleLevel)
{
    string prefix = "inference-";
    string postfix = "CL.txt";
    string filenameString;
    int Ndimensions;
    int Niterations;
    const char *filename = 0;
    ArrayXd marginalParameter;
    ArrayXd marginalDistribution;

    Ndimensions = nestedSampler.posteriorSample.rows();
    Niterations = nestedSampler.posteriorSample.cols();
    inferenceResults.resize(Ndimensions, 5);


    // Now loop over all free parameters

    for (int i = 0; i < Ndimensions; i++)
    {
        marginalParameter = nestedSampler.posteriorSample.row(i);
        marginalDistribution = posteriorDistribution;


        // Sort elements of array marginalParameter in increasing
        // order and sort elements of array marginalDistribution accordingly
        
        Functions::sortElements(marginalParameter, marginalDistribution); 
        

        // Marginalize over parameter by merging posterior values 
        // corresponding to equal parameter values

        int same = 0;

        for (int j = 0; j < Niterations - 1; j++)
        {  
            if (marginalParameter(j) == -DBL_MAX)
                continue;
            else
            {
                for (int k = j + 1; k < Niterations; k++)
                {
                    if (marginalParameter(k) == -DBL_MAX)
                        continue;
                    else
                        if (marginalParameter(j) == marginalParameter(k))
                        {   
                            marginalParameter(k) = -DBL_MAX;        // Set duplicate to bad value (flag)
                            marginalDistribution(j) = marginalDistribution(j) + marginalDistribution(k); // Merge probability values
                            marginalDistribution(k) = 0.0;
                            same++;
                        }
                }
            }
        }


        // Remove bad points and store final values in array copies

        if (same > 0) // Check if bad points are present otherwise skip block
        {
            ArrayXd marginalParameterCopy(Niterations - same);
            ArrayXd marginalDistributionCopy(Niterations - same);
        
            int n = 0;

            for (int m = 0; (m < Niterations) && (n < (Niterations - same)); m++)
            {
                if (marginalParameter(m) == -DBL_MAX)
                    continue;
                else
                    if (marginalParameter(m) != -DBL_MAX)
                        {
                            marginalParameterCopy(n) = marginalParameter(m);
                            marginalDistributionCopy(n) = marginalDistribution(m);
                            n++;
                        }
            }


            // Replace original marginal arrays with array copies

            marginalParameter = marginalParameterCopy;
            marginalDistribution = marginalDistributionCopy;
        }


        /* // BEGIN TEST
        marginalParameter.resize(500);
        marginalDistribution.resize(500);
    
        for (int jj = 0; jj < 500; jj++)
        {
            marginalParameter(jj) = jj/static_cast<double>(marginalParameter.size()) * 30.;
            marginalDistribution(jj) = exp(Functions::logGaussProfile(marginalParameter(jj), 10, 1.5, 1));
        }
        marginalDistribution = marginalDistribution/marginalDistribution.sum();
        // END TEST */ 


        // Compute the mean value (expectation value)
       
        double meanParameter;
        
        meanParameter = (marginalParameter.cwiseProduct(marginalDistribution)).sum();
        inferenceResults(i,0) = meanParameter;


        // Compute the median value (value corresponding to 50% of probability)

        double medianParameter;
        double totalProbability = 0.0;
        int k = 0;
        
        while (totalProbability < 0.5)
        {
            medianParameter = marginalParameter(k);
            totalProbability += marginalDistribution(k);
            k++;
        }
        
        inferenceResults(i,1) = medianParameter;


        // Find the mode value (parameter corresponding to maximum probability value)

        int max = 0;                                    // Subscript corresponding to mode value
        double maximumMarginal;
        double maximumParameter;

        maximumMarginal = marginalDistribution.maxCoeff(&max);
        maximumParameter = marginalParameter(max);
        inferenceResults(i,2) = maximumParameter;

        
        // Compute the "shortest" credible intervals (CI)

        int stepRight = 1;
        double limitProbabilityRight;
        double limitProbabilityLeft;
        double limitParameterLeft;
        double limitParameterRight;
        ArrayXd marginalDistributionLeft(max + 1);       // Marginal distribution up to mode value
        ArrayXd marginalParameterLeft(max + 1);          // Parameter range up to mode value
        ArrayXd marginalDifferenceLeft;                  // Difference distribution to find point in 
                                                         // in the left-hand distribution having equal
                                                         // probability to that identified in the right-hand part
        for (int nn = 0; nn <= max; nn++)
        {
            marginalDistributionLeft(nn) = marginalDistribution(nn);          // Consider left-hand distribution (up to mode value)
            marginalParameterLeft(nn) = marginalParameter(nn);
        }


        while (totalProbability < (credibleLevel/100.))     // Stop when probability >= credibleLevel 
        {
            totalProbability = 0.0;
            limitProbabilityRight = marginalDistribution(max + stepRight);
            limitParameterRight = marginalParameter(max + stepRight);
            marginalDifferenceLeft = (marginalDistributionLeft - limitProbabilityRight).abs();

            int min = 0;

            limitProbabilityLeft = marginalDifferenceLeft.minCoeff(&min);
            limitProbabilityLeft = marginalDistribution(min);
            limitParameterLeft = marginalParameter(min);

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
           
        inferenceResults(i,3) = lowerCredibleInterval;
        inferenceResults(i,4) = upperCredibleInterval;

    }


    // Write output ASCII file

    filenameString = outputDirectory + dataFilename + prefix + to_string(static_cast<int>(credibleLevel)) + postfix;
    filename = filenameString.data();
    ofstream outputFile(filename);

    if (!outputFile.good())
    {
        cerr << "Error opening output file" << endl;
        exit(EXIT_FAILURE);
    }
            
    outputFile << "# Inference results from nested algorithm" << endl;
    outputFile << "# Credible intervals are the shortest credible intervals" << endl; 
    outputFile << "# according to the usual definition" << endl;
    outputFile << "# Credible level: " << fixed << setprecision(2) << credibleLevel << " %" << endl;
    outputFile << "# Column #1: Expectation" << endl;
    outputFile << "# Column #2: Median" << endl;
    outputFile << "# Column #3: Mode" << endl;
    outputFile << "# Column #4: Lower Credible Interval (CI)" << endl;
    outputFile << "# Column #5: Upper Credible Interval (CI)" << endl;
    outputFile << fixed << setprecision(12);
    File::arrayToFile(outputFile, inferenceResults);
    outputFile.close();

} // END Results::printInference()










// Results::getPosteriorDistribution()
//
// PURPOSE:
//      Gets private data member posteriorDistribution.
//
// OUTPUT:
//      An Eigen Array containing the values of the posterior
//      distribution.
// 

ArrayXd Results::getPosteriorDistribution()
{
    return posteriorDistribution;

} // END Results::getPosteriorDistribution()











// Results::getInferenceResults()
//
// PURPOSE:
//      Gets private data member inferenceResults.
//
// OUTPUT:
//      An Eigen Array containing all the values obtained
//      from the inference analysis of the posterior probability
//      distribution.
// 

ArrayXXd Results::getInferenceResults()
{
    return inferenceResults;

} // END Results::getInferenceResults()
