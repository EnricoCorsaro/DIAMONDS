// Class for input/output stream to/from ASCII files
// Enrico Corsaro @ IvS - 22 January 2013
// e-mail: enrico.corsaro@ster,kuleuven.be
// Header file "FileProcess.h"

#ifndef FILEPROCESS_H
#define FILEPROCESS_H

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>


using namespace std;

class FileProcess
{
    
    public:

        unsigned long nlines;
        int ncolumns;
        vector<double> x;
        vector<double> y;
        FileProcess(const char *name);
        void read2ColDouble ();
        
    private:

        unsigned long countLine();
        int countColumn();
        const char *filename;
};



#endif
