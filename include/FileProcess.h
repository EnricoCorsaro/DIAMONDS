// Class for input/output stream to/from ASCII files
// Enrico Corsaro @ IvS - 22 January 2013
// e-mail: enrico.corsaro@ster,kuleuven.be
// Header file "FileProcess.h"

#include "lib.h"
using namespace std;

class FileProcess
{
    public:
        unsigned long nlines;
        int ncolumns;
        vector<double> x;
        vector<double> y;
    
        FileProcess (const char *name)
        {
            filename = name;
            nlines  = countLine();
            ncolumns = countColumn();
        }
    

        // Read an ASCII file in a two columns format containing float numbers
        // and save them in two vectors having double type, both accessible from outside

        void read2ColDouble()
        {
            ifstream inputFile;
            inputFile.open(filename);
        
            if (!inputFile)
            {
                cout << "Could not open file. Quitting program." << endl;
                exit(1);
            }
    
            x.resize(nlines);
            y.resize(nlines);
        
            int i = 0;
            for (int i=0; inputFile && (i < nlines); i++)
            {
                inputFile >> x.at(i) >> y.at(i);
            }
            //inputFile.close();

            return;
        }
    
    
    private:    
    
        // Counts the total number of lines in a file
        
        unsigned long countLine ()  
        {
            ifstream inputFile;
            inputFile.open(filename);

            if (!inputFile)
            {
                cout << "Could not open file. Quitting program." << endl;
                exit(1);
            }

            unsigned long n = 0;
            string s;
            while (getline(inputFile,s)) ++n;

            inputFile.close();
            return n;
        }
    
        // Counts the total number of columns in a file
        
        int countColumn()      
        {
            ifstream inputFile;
            inputFile.open(filename);

            if (!inputFile)
            {
                cout << "Could not open file. Quitting program." << endl;
                exit(1);
            }

            int ncol = 0;
            string line;
            getline(inputFile,line);

            stringstream stringFile(line);
            string ss;
            while (stringFile >> ss) ++ncol;

            inputFile.close();
            return ncol;
        }

        const char *filename;
};
