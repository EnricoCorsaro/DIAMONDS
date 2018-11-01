// Compile with:
// clang++ -o demoPrincipalComponentProjector5D demoPrincipalComponentProjector5D.cpp -L../build/ -I ../include/ -l diamonds -stdlib=libc++ -std=c++11 -Wno-deprecated-register
//

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include "File.h"
#include "PrincipalComponentProjector.h"

using namespace std;
using namespace Eigen;


int main()
{
    // Open the input file and read the data (synthetic sampling of a 5D parameter space)
    
    ifstream inputFile;
    File::openInputFile(inputFile, "kmeans_testsample5D.txt");
    
    unsigned long Nrows;
    int Ncols;

    File::sniffFile(inputFile, Nrows, Ncols);
    ArrayXXd data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    ArrayXXd sample = data.transpose();
    inputFile.close();

    
    // Set up the Principal Component Projector and apply this to the input sample

    bool printNdimensions = true;
    PrincipalComponentProjector projector(printNdimensions);
    
    ArrayXXd optimizedSample;
    optimizedSample = projector.projection(sample);
    ArrayXXd finalSample = optimizedSample.transpose();


    // Print the reduced-dimensionality sample into an output ASCII file

    ofstream outputFile;
    File::openOutputFile(outputFile, "principalComponentProjection5D.txt");
    outputFile << scientific << setprecision(4);
    File::arrayXXdToFile(outputFile, finalSample);
    outputFile.close();



    // That's it!
 
    return EXIT_SUCCESS;
}
