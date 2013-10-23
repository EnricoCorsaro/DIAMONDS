// Namespace for file processing.
// Joris De Ridder & Enrico Corsaro @ IvS - 12 February 2013
// e-mail: joris.deridder@ster.kuleuven.be
// Header file "File.h"
// Implementation contained in "File.cpp"


#ifndef FILE_H
#define FILE_H


#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <Eigen/Core>


using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


namespace File
{
    void openInputFile(ifstream &inputFile, string inputFileName);
    void openOutputFile(ofstream &outputFile, string outputFileName);
    
    ArrayXXd arrayXXdFromFile(ifstream &inputFile, const unsigned long Nrows, const int Ncols, char separator = ' ', char commentChar = '#');
    void arrayXXdToFile(ofstream &outputFile, RefArrayXXd array, string separator = "  ", string terminator = "\n");
    void twoArrayXdToFile(ofstream &outputFile, RefArrayXd array1, RefArrayXd array2, string separator = "  ", string terminator = "\n");
    void arrayXdToFile(ofstream &outputFile, RefArrayXd array, string terminator = "\n");
    void arrayXXdRowsToFiles(RefArrayXXd array, string fullPathPrefix, string fileExtension = ".txt", string terminator = "\n");
    void sniffFile(ifstream &inputFile, unsigned long &Nrows, int &Ncols, char separator = ' ', char commentChar = '#');
    void configuringParametersToFile(ofstream &outputFile, const int initialNobjects, const int minNobjects, 
                                     const int minNclusters, const int maxNclusters, const int NinitialIterationsWithoutClustering,
                                     const int NiterationsWithSameClustering, const int maxNdrawAttempts, 
                                     const double initialEnlargementFraction, const double shrinkingRate, const double terminationFactor);

} // END namespace File

#endif
