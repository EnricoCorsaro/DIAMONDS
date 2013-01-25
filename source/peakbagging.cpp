// Main code for peak bagging by means of nested sampling analysis
// Created by Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Source code file "peakbagging.cpp"

#include <iostream>
#include <iomanip>
#include "MathExtra.h"
#include "FileProcess.h"
#include "NestedSampler.h"

int main()
{
    int Ndata = 100;      // Number of objects per nested iteration (usually 100)
    int Niter = 1000;     // Number of nested iterations (usually 1000)
    
    NestedSampler nestedSampler(Ndata, Niter);
    nestedSampler.run();
    
    cout << right << setw(10) << "Parameter value" << right << setw(20) << "logLikelihood" << endl;
    for (int i = 0; i < Niter; i++)
    {
       cout << right << setw(10) << nestedSampler.postP[i]
       << right << setw(20) << nestedSampler.postlogL[i] << endl;
    }
    
    cout << endl;
    cout << " ------------------------------------------------" << endl;
    cout << " Evidence: logZ = " << nestedSampler.results[0] << " +/- " << nestedSampler.results[1] << endl;
    cout << " Information: H = " << nestedSampler.results[2] << endl;
    cout << " ------------------------------------------------" << endl;
}
