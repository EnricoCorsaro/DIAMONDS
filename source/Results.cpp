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
    string name = "loglikelihood.txt";
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









// Results::printPosterior()
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

void Results::printPosterior()
{
    ArrayXd logPosterior;
    string name = "posterior.txt";
    string filenameString;
    const char *filename = 0;

    logPosterior = nestedSampler.logWeightOfPosteriorSample 
    - nestedSampler.getLogEvidence();
    posteriorOfPosteriorSample = logPosterior.exp();

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
    File::oneArrayToFile(outputFile, posteriorOfPosteriorSample);
    outputFile.close();

} // END Results::printPosterior()










// Results::printInference()
//
// PURPOSE:
//      Prints the expectation values from the marginalized posterior 
//      probability density into an ASCII file. Bayesian credible 
//      intervals (CI) are also computed.
//      The values are also stored into the one dimensional Eigen Array
//      PosteriorOfPosteriorSample.
//
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
    string name = "inference.txt";
    string filenameString;
    int Ndimensions;
    double meanValue;
    const char *filename = 0;
    double posteriorParameterMaximum;
    ArrayXd posteriorParameter;
    ArrayXd weightedValues;

    Ndimensions = nestedSampler.posteriorSample.rows();
    ArrayXXd expectations(Ndimensions, 3);

    for (int i = 0; i < Ndimensions; i++)
    {
        // Compute the mean value (expectation value)
       
        posteriorParameter = nestedSampler.posteriorSample.row(i);
        meanValue = (posteriorParameter.cwiseProduct(posteriorOfPosteriorSample)).sum();
        expectations(i,0) = meanValue;


        // Find element corresponding to maximum posterior

        int max = 0;
        double maximumPosterior;

        maximumPosterior = posteriorOfPosteriorSample.maxCoeff(&max);
        posteriorParameterMaximum = posteriorParameter(max);

        /*
        // Compute the credible intervals
        // Now sorting posteriorParameter into increasing order
        // and sorting posteriorOfPosteriorSample according to new
        // sorting of posteriorParameter elements

        ArrayXd posteriorParameterSorted = posteriorParameter;
        ArrayXd posteriorOfPosteriorSampleSorted = posteriorOfPosteriorSample;
        
        MathExtra::sortElements(posteriorParameterSorted, posteriorOfPosteriorSampleSorted); 
        double probability = 0.0;   // Start with zero probability
        int step = 0;

        do
        {
            
            ++step;    
        }
        while (probability <= 0.01*credibleLevel); // Stop when probability is <= credibleLevel */
    }
    filenameString = outputDirectory + name;
    filename = filenameString.data();
    ofstream outputFile(filename);

    if (!outputFile.good())
    {
        cerr << "Error opening output file" << endl;
        exit(EXIT_FAILURE);
    }
            
    outputFile << "# Expectation values from nested algorithm" << endl;
    outputFile << "# Credible level: " << fixed << setprecision(2) << credibleLevel << " %" << endl;
    outputFile << "# Column #1: Expectations" << endl;
    outputFile << "# Column #2: Lower Credible Intervals (CI)" << endl;
    outputFile << "# Column #3: Upper Credible Intervals (CI)" << endl;
    outputFile << fixed << setprecision(12);
    File::arrayToFile(outputFile, expectations);
    outputFile.close();

} // END Results::printInference()
