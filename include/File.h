// Namespace for file processing.
// Joris De Ridder & Enrico Corsaro @ IvS - 12 February 2013
// e-mail: joris.deridder@ster.kuleuven.be / emncorsaro@gmail.com
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
    vector<string> vectorStringFromFile(ifstream &inputFile, const unsigned long Nrows, char commentChar = '#');

    void arrayXXdToFile(ofstream &outputFile, RefArrayXXd array, string separator = "  ", string terminator = "\n");
    void twoArrayXdToFile(ofstream &outputFile, RefArrayXd array1, RefArrayXd array2, string separator = "  ", string terminator = "\n");
    void arrayXdToFile(ofstream &outputFile, RefArrayXd array, string terminator = "\n");
    void arrayXXdRowsToFiles(RefArrayXXd array, string fullPathPrefix, string fileExtension = ".txt", string terminator = "\n");
    void sniffFile(ifstream &inputFile, unsigned long &Nrows, int &Ncols, char separator = ' ', char commentChar = '#');

}

#endif
