
#include "File.h"



// File::openInputFile()
//
// PURPOSE: 
//      Opens a file to read from and check its sanity
//
// INPUT:
//      inputFile: clean input stream
//      inputFileName: full path of the input file.
// 
// OUTPUT:
//      void
//

void File::openInputFile(ifstream &inputFile, string inputFileName)
{
    inputFile.open(inputFileName.c_str());
    if (!inputFile.good())
    {
        cerr << "Error opening input file " << inputFileName << endl;
        exit(EXIT_FAILURE);
    }
}









// File::openOutputFile()
//
// PURPOSE: 
//      Opens a file to write and check its sanity
//
// INPUT:
//      outputFile: clean output stream
//      outputFileName: full path of the output file.
// 
// OUTPUT:
//      void
//

void File::openOutputFile(ofstream &outputFile, string outputFileName)
{
    outputFile.open(outputFileName.c_str());
    if (!outputFile.good())
    {
        cerr << "Error opening input file " << outputFileName << endl;
        exit(EXIT_FAILURE);
    }
}











// File::arrayXXdFromFile()
//
// PURPOSE: 
//      Reads an ascii file into an Eigen ArrayXXd
//
// INPUT:
//      inputFile: input stream. Supposed to be already opened, and checked for sanity.
//      Nrows: number of rows with values in the file
//      Ncols: number of columns in each row
//      separator: numbers are normally separated by ' '. Use ',' for CSV files.
//      commentChar: all lines starting with this character (e.g. '#') are skipped
// 
// OUTPUT:
//      An Eigen::ArrayXXd.
//
//  REMARKS:
//      - the stream is not closed afterwards.
//

ArrayXXd File::arrayXXdFromFile(ifstream &inputFile, const unsigned long Nrows, const int Ncols, char separator, char commentChar)
{
    string line;
    unsigned long iRow = 0;
    int iCol = 0;
    ArrayXXd array(Nrows, Ncols);
  
    while(!inputFile.eof())
    {
        getline(inputFile, line);
    

        // Skip those lines that start with the comment character
        
        if (line[0] == commentChar) continue; 
    

        // Skip those lines with only whitespace
        
        const std::string whitespace = " \t";
        const auto strBegin = line.find_first_not_of(whitespace);
        if (strBegin == string::npos) continue;
    

        // Check if the line number doesn't exceed the expected number of lines
        
        if (iRow > Nrows-1)
        {
            cerr << "Error: numbers of rows in file exceeds " << Nrows << endl;
            exit(EXIT_FAILURE);
        }
    

        // Tokenize the line into chuncks delimited by separator
        
        vector<string> tokens;
        string::size_type begin = line.find_first_not_of(separator);
        string::size_type end = line.find_first_of(separator, begin);
        while (begin != string::npos || end != string::npos)
        {
            tokens.push_back(line.substr(begin, end-begin));
            begin = line.find_first_not_of(separator, end);
            end = line.find_first_of(separator, begin);
        }
    

        // Check if the number of numbers on the line matches the expected
        // number of columns
        
        if (tokens.size() != Ncols)
        {
            cerr << "Error on row " << iRow << ": number of tokens != " 
                 << Ncols << endl;
            exit(EXIT_FAILURE);
        }
    

        // Convert the strings to numbers, and store them in the array
        
        iCol = 0;
        for (auto token : tokens)
        {
            istringstream mystream(token);
            double value;
            if (!(mystream >> value))
            {
                cerr << "Error on row " << iRow << ". Can't convert "
                     << token << " to a number" << endl;
                exit(EXIT_FAILURE);
            }
            else
            {
                array(iRow, iCol) = value;
                iCol++;
            }
        }
        
        iRow++;
    }
    

    // Check if the current line number is less than the expected number of lines
    
    if (iRow < Nrows)
    {
        cerr << "Warning: too few lines in the file: " << iRow << " < " << Nrows << endl;
    }
    
    return array;
}


    








// File::arrayXXdToFile()
//
// PURPOSE: writes an Eigen ArrayXXd to an ascii file
//
// INPUT:
//      outputFile: output stream, assumed to be already opened and checked for sanity.
//      array: array to write to a file
//      separator: column separator, e.g. "  "
//      terminator: line terminator, e.g. "\n"
// 
// OUTPUT:
//      - file contents is changed
// 
// REMARKS:
//      - the stream is not closed afterwards
//      - this function can equally well be used to append to an existing file
//

void File::arrayXXdToFile(ofstream &outputFile, ArrayXXd array, string separator, string terminator)
{
    for (ptrdiff_t i = 0; i < array.rows(); ++i)
    {
        for (ptrdiff_t j = 0; j < array.cols()-1; ++j)
        {
            outputFile << array(i,j) << separator;
        }
        outputFile << array(i,array.cols()-1) << terminator;
    }
} 








// File::twoArrayXdToFile()
//
// PURPOSE: 
//      Writes two Eigen::ArrayXd arrays as columns to an ascii file
//
// INPUT:
//      outputFile: output stream, assumed to be already opened and checked for sanity.
//      array1: array which will appear as the first column in the file
//      array2: array which will appear as the first column in the file
//      separator: column separator, e.g. "  "
//      terminator: line terminator, e.g. "\n"
// 
// OUTPUT:
//      - file contents is changed
// 
// REMARKS:
//      - overloaded function
//      - the stream is not closed afterwards
//      - this function can equally well be used to append to an existing file
//

void File::twoArrayXdToFile(ofstream &outputFile, ArrayXd array1, ArrayXd array2, string separator, string terminator)
{
    assert(array1.size() == array2.size());
    
    for (ptrdiff_t i = 0; i < array1.rows(); ++i)
    {
        outputFile << array1(i) << separator << array2(i) << terminator;
    }
}











// File::arrayXdToFile()
//
// PURPOSE: 
//      Writes onw Eigen::ArrayXd arrays as a column to an ascii file
//
// INPUT:
//      outputFile: output stream, assumed to be already opened and checked for sanity.
//      array: array which will appear as a column in the file
//      terminator: line terminator, e.g. "\n"
// 
// OUTPUT:
//      - file contents is changed
// 
// REMARKS:
//      - the stream is not closed afterwards
//      - this function can equally well be used to append to an existing file
//

void File::arrayXdToFile(ofstream &outputFile, ArrayXd array, string terminator)
{
    for (ptrdiff_t i = 0; i < array.size(); ++i)
    {
        outputFile << array(i) << terminator;
    }
} 











void File::arrayXXdRowsToFiles(ArrayXXd array, string fullPathPrefix, string fileExtension, string terminator)
{
    int Nrows = array.rows();
    assert(Nrows > 0);


    // Find out the number of decimal digits that the number of dimensions has
    
    int Ndigits = int(floor(log10(double(Nrows)))); 
    

    // Write everything to the output files

    for (int i = 0; i < Nrows; i++)
    {
        // Include the row number with preceding zeros in the filename
        
        ostringstream numberString;
        numberString << setfill('0') << setw(3-Ndigits) << i;
        string fullPath = fullPathPrefix + numberString.str() + fileExtension;
       

        // Open the output file and check for sanity
        
        ofstream outputFile;
        File::openOutputFile(outputFile, fullPath);
       

        // Write all values of this particular parameter in our sample to the output file
        
        outputFile << setiosflags(ios::scientific) << setprecision(9);
        File::arrayXdToFile(outputFile, array.row(i), terminator);
        outputFile.close();
    }
}











// File::snifFile()
//
// PURPOSE: 
//      Checks a file and derives the number of rows and columns
//
// INPUT:
//      inputFile: input stream, assumed to be already opened and checked for sanity.
//      separator: column separator, e.g. "  "
//      Nrows: see output
//      Ncols: see output
// 
// OUTPUT:
//      - Nrows will contain the number of rows
//      - Ncols will constina the number of columns
//
// REMARKS:
//      - the stream is not closed afterwards, but it is rewinded.
//

void File::snifFile(ifstream &inputFile, unsigned long &Nrows, int &Ncols, char separator, char commentChar)
{
    string line;
    int iRow = 0;
    
    while(!inputFile.eof())
    {
        getline(inputFile, line);
    

        // Skip those lines that start with the comment character
        
        if (line[0] == commentChar) continue; 
    

        // Skip those lines with only whitespace
        
        const std::string whitespace = " \t";
        const auto strBegin = line.find_first_not_of(whitespace);
        if (strBegin == string::npos) continue;
        
                
        // For the very first line, tokenize the line into chuncks delimited by separator
        
        if (iRow == 0)
        {
            vector<string> tokens;
            string::size_type begin = line.find_first_not_of(separator);
            string::size_type end = line.find_first_of(separator, begin);
            while (begin != string::npos || end != string::npos)
            {
                tokens.push_back(line.substr(begin, end-begin));
                begin = line.find_first_not_of(separator, end);
                end = line.find_first_of(separator, begin);
            }
            
            Ncols = tokens.size();
        }
                
        iRow++;
    }
    
    Nrows = iRow;
    

    // Rewind the input stream to the beginning
    
    inputFile.clear();
    inputFile.seekg(ios::beg);
}



