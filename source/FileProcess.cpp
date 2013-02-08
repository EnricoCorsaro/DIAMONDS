#include "FileProcess.h"


// FileProcess::FileProcess()
// 
// PURPOSE: 
//      Constructor
// INPUT:
//      Filename to read/write
// OUTPUT:
//

FileProcess::FileProcess(const char *name)
{
    filename = name;
}





// FileProcess::read2ColDouble()
// PURPOSE: 
//      Reads an ASCII file having a two columns format of doubles and stores each
//      column in a vector returned as a public data member
// INPUT:
// OUTPUT:
//

void FileProcess::read2ColDouble()
{
    nlines = countLine();
    ncolumns = countColumn();
    
    if (ncolumns > 2)
    {
        cout << "Number of columns exceeds the expectations. Quitting program." << endl;
        exit(1);
    }
    
    ifstream inputFile;
    inputFile.open(filename);

    if (!inputFile)
    {
        cout << "Could not open file. Quitting program." << endl;
        exit(1);
    }

    x.resize(nlines);
    y.resize(nlines);

    for (unsigned long i = 0; i < nlines; i++)
    {
        if (inputFile.good() && !inputFile.eof())
        {
            inputFile >> x.at(i) >> y.at(i);
        }
    }
    inputFile.close();

    return;
}




// FileProcess::write2ColDouble()
// PURPOSE: 
//      Writes an ASCII file having a two columns format of doubles taken from
//      two input vectors.
// INPUT:
//      x1: vector of independent values
//      x2: vector of dependent values
//      title1: string containing label for x1 (default blank)
//      title2: string containing label for x2 (default blank)
// OUTPUT:
//

void FileProcess::write2ColDouble(vector<double> &x1, vector<double> &x2, string title1, string title2)
{
    unsigned long x1size,x2size;
    x1size = x1.size();
    x2size = x2.size();

    if (x1size != x2size)
    {
        cout << "Array dimensions do not match. Quitting program." << endl;
        exit(1);
    }

    ofstream outputFile;
    outputFile.open(filename);

    if (!outputFile)
    {
        cout << "Could not open file. Quitting program." << endl;
        exit(1);
    }

    outputFile << left << "#" << right << setw(12) << title1 << right << setw(22) << title2 << endl;

    for (unsigned long i = 0; i < x1size; i++)
    {
        if (outputFile.good() && !outputFile.eof())
        {
            outputFile << right << setw(12) << setprecision(8) << x1[i]
            << right << setw(22) << setprecision(8) << x2[i] 
            << endl;
        }
    }
    outputFile.close();

    return;
}




// FileProcess::countLine()
// PURPOSE: 
//      Counts the total number of lines in a file
// INPUT:
// OUTPUT:
//      An unsigned long containing the number of lines.
//

unsigned long FileProcess::countLine()
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

    while (getline(inputFile,s))
        ++n;

    inputFile.close();
    return n;
}





// FileProcess::countColumn()
// PURPOSE: 
//      Counts the total number of columns in a file
// INPUT:
// OUTPUT:
//      An integer containing the number of columns
//

int FileProcess::countColumn ()
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

    while (stringFile >> ss)
        ++ncol;

    inputFile.close();
    return ncol;
}
