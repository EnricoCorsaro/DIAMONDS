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
#include <fstream>
#include <sstream>
#include <Eigen/Core>

using namespace std;
using Eigen::ArrayXXd;
using Eigen::ArrayXd;

namespace File
{
    ArrayXXd arrayFromFile(ifstream &inputFile, const unsigned long Nrows, const int Ncols, char separator = ' ', char commentChar = '#');
    void arrayToFile(ofstream &outputFile, ArrayXXd array, string separator = "  ", string terminator = "\n");
    void arrayToFile(ofstream &outputFile, ArrayXd array1, ArrayXd array2, string separator = "  ", string terminator = "\n");
    void OneArrayToFile(ofstream &outputFile, ArrayXd array, string terminator = "\n");
    void snifFile(ifstream &inputFile, unsigned long &Nrows, int &Ncols, char separator = ' ', char commentChar = '#');
} // END namespace File

#endif
