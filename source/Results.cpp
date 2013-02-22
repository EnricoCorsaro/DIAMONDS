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

Results::Results(const NestedSampler &nestedSampler, const char *outputDirectory)
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
            File::OneArrayToFile(outputFile, nestedSampler.posteriorSample.row(i));
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
                File::OneArrayToFile(outputFile, nestedSampler.posteriorSample.row(i));
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
                File::OneArrayToFile(outputFile, nestedSampler.posteriorSample.row(i));
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
    string name = "loglikelihood";
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
    File::OneArrayToFile(outputFile, nestedSampler.logLikelihoodOfPosteriorSample);
    outputFile.close();
} // END Results::printLogLikelihood()
