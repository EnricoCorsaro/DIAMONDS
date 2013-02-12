
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
    ArrayXXd arrayFromFile(ifstream &inputFile, int Nrows, int Ncols, char separator = ' ', char commentChar = '#');
    void arrayToFile(ofstream &outputFile, ArrayXXd array, string separator = "  ", string terminator = "\n");
    void arrayToFile(ofstream &outputFile, ArrayXd array1, ArrayXd array2, string separator = "  ", string terminator = "\n");
    void snifFile(ifstream &inputFile, int &Nrows, int &Ncols, char separator = ' ', char commentChar = '#');
}


#endif
