// Demonstration of classes included as header files
// Created by Enrico Corsaro @ IvS - 22 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be

#include "../include/MathExtra.h"
#include "../include/FileProcess.h"

int main()
{
    vector<double> x(100);
    vector<double> y(100);
    double amp, gamma, x0;

    MathExtra math1;

    for (int i=0; i < x.size(); i++)
    {
        x.at(i) = i;
        y.at(i) = 2.;
    }

    amp = 2.;
    gamma = 3.;
    x0 = 20.;

    cout << math1.total(x) << endl;
    cout << setprecision(20) << math1.product(y) << endl;

    //math1.lorentzProfile(amp,gamma,x0,x,y)    ;


//  Write output file
    
//  ofstream outputFile         ;
//  outputFile.open ("test.dat", ofstream::out) ;   

//  if (!outputFile)
//  {
//      cout << "Could not open file. Closing session." << endl;
//      exit(1);
//  }
    
//  outputFile << "x" << setw(10) << "y" << endl;

//  for (int i=0; i < x.size(); i++)
//  {
//      outputFile << setprecision(2) << fixed << x.at(i) << setw(10) << setprecision(6) 
//          << fixed << y.at(i) << endl;
//  }
//  outputFile.close();
    
    double  xin, yin;
    unsigned long nlines;
    int ncolumns;
    //const char * filename = "test.dat"        ;
    //FileProcess   file1 ( "test.dat" )        ;

    return 0;
}
