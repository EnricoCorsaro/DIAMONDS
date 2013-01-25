// Main code for peak bagging by means of nested sampling analysis
// Created by Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Source code file "peakbagging.cpp"
#include "../include/Nesting.h"

int main()
{
    int n_obj;      // Number of objects per nested iteration (usually 100)
    int n_nest;     // Number of nested iterations (usually 1000)
    
    n_obj = 100;
    n_nest = 1000;
    Nesting nestObj( n_obj, n_nest );
    
    cout << right << setw(10) << "Parameter value" << right << setw(20) << "logLikelihood" << endl;
    for (int i = 0; i < n_nest; i++ )
    {
       cout << right << setw(10) << nestObj.postP.at(i) 
       << right << setw(20) << nestObj.postlogL.at(i) << endl;
    }
    
    cout << endl;
    cout << " ------------------------------------------------" << endl;
    cout << " Evidence: logZ = " << nestObj.results.at(0) << " +/- " << nestObj.results.at(1) << endl;
    cout << " Information: H = " << nestObj.results.at(2) << endl;
    cout << " ------------------------------------------------" << endl;
}
