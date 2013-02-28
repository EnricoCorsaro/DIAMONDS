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

Results::Results(NestedSampler &nestedSampler, const char *outputDirectory)
: outputDirectory(outputDirectory),
  nestedSampler(nestedSampler)
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
            filenameString = outputDirectory + prefix1 + to_string(i) + postfix;
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
                filenameString = outputDirectory + prefix2 + to_string(i) + postfix;
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
                filenameString = outputDirectory + prefix3 + to_string(i) + postfix;
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
    string name = "log-likelihood.txt";
    string filenameString;
    const char *filename = 0;

    filenameString = outputDirectory + name;
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

    filenameString = outputDirectory + name;
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









// Results::printPosteriorDensity()
//
// PURPOSE:
//      Prints the posterior probability density for the sample 
//      obtained from the nested sampling into an ASCII file of one column format. 
//      The values are also stored into the one dimensional Eigen Array
//      PosteriorOfPosteriorSample.
//
// OUTPUT:
//      void
// 

void Results::printPosteriorDensity()
{
    ArrayXd logPosteriorDensity;
    string name = "posterior-density.txt";
    string filenameString;
    const char *filename = 0;

    logPosteriorDensity = nestedSampler.logWeightOfPosteriorSample 
    - nestedSampler.getLogEvidence();
    posteriorDensity = logPosteriorDensity.exp();

    filenameString = outputDirectory + name;
    filename = filenameString.data();
    ofstream outputFile(filename);

    if (!outputFile.good())
    {
        cerr << "Error opening output file" << endl;
        exit(EXIT_FAILURE);
    }
            
    outputFile << "# Posterior probability density from nested algorithm" << endl;
    outputFile << scientific << setprecision(12);
    File::oneArrayToFile(outputFile, posteriorDensity);
    outputFile.close();

} // END Results::printPosteriorDensity()










// Results::printInference()
//
// PURPOSE:
//      Prints the expectation, median and mode values from the 
//      marginalized posterior probability density into an ASCII file. 
//      Shortest Bayesian credible intervals (CI) are also computed.
//      All the values are stored in to the Eigen Array inference.

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
    string prefix = "inference_";
    string postfix = "CL.txt";
    string filenameString;
    int Ndimensions;
    int Niterations;
    double meanValue;
    const char *filename = 0;
    ArrayXd marginalParameter;
    ArrayXd marginalDensity;

    Ndimensions = nestedSampler.posteriorSample.rows();
    Niterations = nestedSampler.posteriorSample.cols();
    inference.resize(Ndimensions, 5);
    


    // Now loop over all free parameters

    for (int i = 0; i < Ndimensions; i++)
    {
        marginalParameter = nestedSampler.posteriorSample.row(i);
        marginalDensity = posteriorDensity;


        // Sort elements of array marginalParameter in increasing
        // order and those of array marginalDensity accordingly
        
        MathExtra::sortElements(marginalParameter, marginalDensity); 
        

        // Marginalize over parameter by merging posterior densities of 
        // equal parameter values

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
                            marginalDensity(j) = marginalDensity(j)+marginalDensity(k); // Merge probability densities
                            marginalDensity(k) = 0.0;
                            same++;
                        }
                }
            }
        }


        // Remove bad points and store final values in array copies

        ArrayXd marginalParameterCopy(Niterations - same);
        ArrayXd marginalDensityCopy(Niterations - same);
        
        int n = 0;

        for (int m = 0; (m < Niterations) && (n < (Niterations - same)); m++)
        {
            if (marginalParameter(m) == -DBL_MAX)
                continue;
            else
                if (marginalParameter(m) != - DBL_MAX)
                    {
                        marginalParameterCopy(n) = marginalParameter(m);
                        marginalDensityCopy(n) = marginalDensity(m);
                        n++;
                    }
        }


        // Replace original marginal arrays with array copies

        marginalParameter = marginalParameterCopy;
        marginalDensity = marginalDensityCopy;
        

        // Compute the mean value (expectation value)
       
        meanValue = (marginalParameter.cwiseProduct(marginalDensity)).sum();
        inference(i,0) = meanValue;


        // Compute the median value (value corresponding to 50% of probability)

        double medianValue;
        double totalProbability = 0.0;
        int k = 0;
        
        while (totalProbability < 0.5)
        {
            medianValue = marginalParameter(k);
            totalProbability += marginalDensity(k);
            k++;
        }
        
        inference(i,1) = medianValue;


        // Find the mode value (parameter corresponding to maximum probability density value)

        int max = 0;
        double maximumMarginal;

        maximumMarginal = marginalDensity.maxCoeff(&max);
        inference(i,2) = marginalParameter(max);

        
        /*// Compute the credible intervals

        int step = 0;

        totalProbability = 0.0;

        do
        {
            
            ++step;    
        }
        while (probability < 0.01*credibleLevel); // Stop when probability is >= credibleLevel */
    }


    // Write output ASCII file

    filenameString = outputDirectory + prefix + to_string(static_cast<int>(credibleLevel)) + postfix;
    filename = filenameString.data();
    ofstream outputFile(filename);

    if (!outputFile.good())
    {
        cerr << "Error opening output file" << endl;
        exit(EXIT_FAILURE);
    }
            
    outputFile << "# Inference results from nested algorithm" << endl;
    outputFile << "# Credible intervals are shortest credible intervals according to usual definition" << endl;
    outputFile << "# Credible level: " << fixed << setprecision(2) << credibleLevel << " %" << endl;
    outputFile << "# Column #1: Expectation" << endl;
    outputFile << "# Column #2: Median" << endl;
    outputFile << "# Column #3: Mode" << endl;
    outputFile << "# Column #4: Lower Credible Interval (CI)" << endl;
    outputFile << "# Column #5: Upper Credible Interval (CI)" << endl;
    outputFile << fixed << setprecision(12);
    File::arrayToFile(outputFile, inference);
    outputFile.close();

} // END Results::printInference()










// Results::getPosteriorDensity()
//
// PURPOSE:
//      Gets private data member posteriorDensity.
//
// OUTPUT:
//      An Eigen Array containing the values of the posterior
//      density distribution.
// 

ArrayXd Results::getPosteriorDensity()
{
    return posteriorDensity;

} // END Results::getPosteriorDensity()











// Results::getInference()
//
// PURPOSE:
//      Gets private data member inference.
//
// OUTPUT:
//      An Eigen Array containing all the values obtained
//      from the inference analysis of the posterior density.
// 

ArrayXXd Results::getInference()
{
    return inference;

} // END Results::getInference()
