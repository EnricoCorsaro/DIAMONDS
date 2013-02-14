// Main code for peak bagging by means of nested sampling analysis
// Created by Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Source code file "peakbagging.cpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "MathExtra.h"
#include "File.h"
#include "NestedSampler.h"

int main(int argc, char *argv[])
{
    unsigned long Nrows;
    int Ncols
    ArrayXXd data;
    
    if (argc != 2)
    {
        cerr << "Usage: peakbagging <inputFile> <outputFile>" << endl;
        exit(EXIT_FAILURE);
    }

    ifstream inputFile(argv[1]);
    if (!inputFile.good())
    {
        cerr << "Error opening input file" << endl;
        exit(EXIT_FAILURE);
    }

    File::snifFile(inputFile, Nrows, Ncols);
    data = File::arrayFromFile(inputFile, Nrows, Ncols);

    ofstream outputFile(argv[2]);
    if (!outputFile.good())
    {
        cerr << "Error opening output file" << endl;
        exit(EXIT_FAILURE);
    }
       
    int Nobjects = 100;      // Number of objects per nested iteration (usually 100)
    int Niter    = 1000;     // Number of nested iterations (usually 1000)
    int Ndim     = 1;        // Number of free parameters (dimensions) of the problem
    Prior prior(Ndim);
    ArrayXd parametersMin(Ndim);
    ArrayXd parameterspMax(Ndim);

    // Should give input values from file
    parametersMin = 0.0;
    parametersMax = 20.0;

    prior.setBoundaries(parametersMin, parametersMax);
    
    NestedSampler nestedSampler(Ndim);
    nestedSampler.run(Nobjects, Niter);

    outputFile << "# Parameter value    logLikelihood" << endl;
    outputFile << setiosflags(ios::fixed) << setprecision(8);
    File::arrayToFile(outputFile, nestedSampler.posteriorSample.row(0), nestedSampler.logLikelihoodOfPosteriorSample);
    outputFile.close();
    
    cerr << " Evidence: logZ = " << nestedSampler.getLogEvidence() << " +/- " << nestedSampler.getLogEvidenceError() << endl;
    cerr << " Information: H = " << nestedSampler.getInformationH() << endl;
    
    return EXIT_SUCCESS;
}
