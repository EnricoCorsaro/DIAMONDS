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
#include "NormalVariate.h"


int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cerr << "Usage: peakbagging <outputfile>" << endl;
        exit(EXIT_FAILURE);
    }

    ofstream outputFile(argv[1]);
    if (!outputFile.good())
    {
        cerr << "Error opening output file" << endl;
        exit(EXIT_FAILURE);
    }
       

    int Nobjects = 100;      // Number of objects per nested iteration (usually 100)
    int Niter    = 1000;     // Number of nested iterations (usually 1000)
    
    NormalVariate normalVariate(10.0, 3.0);
    normalVariate.setBoundaries(0.0, 20.0);
    NestedSampler nestedSampler(normalVariate);
    nestedSampler.run(Nobjects, Niter);

    outputFile << "# Parameter value    logLikelihood" << endl;
    outputFile << setiosflags(ios::fixed) << setprecision(8);
    File::arrayToFile(outputFile, nestedSampler.posteriorSample.row(0), nestedSampler.logLikelihoodOfPosteriorSample);
    outputFile.close();
    
    cerr << " Evidence: logZ = " << nestedSampler.getLogEvidence() << " +/- " << nestedSampler.getLogEvidenceError() << endl;
    cerr << " Information: H = " << nestedSampler.getInformationH() << endl;
    
    return EXIT_SUCCESS;
}