
#include "FileProcess.h"


// FileProcess::FileProcess()
// 
// PURPOSE: constructor
//
// INPUT:
//
// OUTPUT:
//

FileProcess::FileProcess(const char *name)
{
    filename = name;
    nlines = countLine();
    ncolumns = countColumn();
}





// FileProcess::read2ColDouble()
// 
// PURPOSE: Reads an ASCII file having a two columns format of doubles and stores each
//          column in a vector returned as a public data member
//
// INPUT:
//
// OUTPUT:
//

void FileProcess::read2ColDouble()
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

    for (unsigned int i = 0; i < nlines; i++)
    {
        if (inputFile.good() && !inputFile.eof())
        {
            inputFile >> x.at(i) >> y.at(i);
        }
    }
    inputFile.close();

    return;
}






// FileProcess::countLine()
// 
// PURPOSE: Counts the total number of lines in a file
//
// INPUT:
//
// OUTPUT:
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
// 
// PURPOSE: Counts the total number of columns in a file
//
// INPUT:
//
// OUTPUT:
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
