#include "MixedPriorMaker.h"


// Functions::preparePriorDistributions()
// 
// PURPOSE:
//      This function builds up mixed prior distributions by reading the hyper parameters and prior types from
//      an input file. The file has to be in a five-columns format, where the first four are related to the prior 
//      hyper parameters and the last one to the prior type. 
//      The prior type is specified through a flag: 0 for Uniform, 1 for Gaussian, 2 for Super Gaussian, 
//      3 for Grid Uniform. Each line of the file will refer to a different free parameter of the model and
//      it can be specified with its own prior type. The prior hyper parameters are defined in the order listed below.
//
//      Uniform prior:
//      Col 1: lower boundary, Col 2: upper boundary, Col 3: 0, Col 4: 0
//
//      Gaussian prior:
//      Col 1: mean value, Col 2: standard deviation, Col 3: 0, Col 4: 0
//
//      Super Gaussian prior:
//      Col 1: mean value, Col 2: standard deviation, Col 3: width of plateau
//
//      Grid Uniform prior: 
//      Col 1: starting point of the grid, Col 2: ending point of the grid, Col 3: number of grid points, Col 4: tolerance (0-1)
//
// INPUT:
//      inputFileName:                   a string containing the name (with full path) for the file containing the prior
//                                       hyper parameters and types
//      outputPathPrefix:                a string containing the full path of the output file of the hyper parameters
//      Ndimensions:                     an empty integer variable to contain the number of free parameters identified
//      writePriorHyperParametersToFile: a boolean to specify whether or not the prior hyper parameters (and prior types)
//                                       must be stored in an output ASCII file
//
// OUTPUT:
//      A vector of pointers to prior class object that will contain the prior distributions 
//
// REMARK:
//      The function operates also with a two-column format input file in the case of uniform prior distribution 
//      for all the free parameters.
//      

vector<Prior*> MixedPriorMaker::prepareDistributions(string inputFileName, string outputPathPrefix, unsigned long &Ndimensions, 
                                                          bool writePriorHyperParametersToFile)
{
    // Read prior hyper parameters for resolved modes
    int Ncols;
    string fullPathOutputHyperParameters;
    ifstream inputFile;
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Ndimensions, Ncols);

    
    ArrayXXd hyperParameters;
    hyperParameters = File::arrayXXdFromFile(inputFile, Ndimensions, Ncols);
    inputFile.close(); 

    if (Ncols != 2)
    {
        if (Ncols == 5)
        {
            cout << "------------------------------------------------------- " << endl;
            cout << " Reading mixed prior types file " + inputFileName << endl;
            cout << "------------------------------------------------------- " << endl;
            cout << endl;

            // Count the number of prior blocks (where each block contains only one prior type) and the corresponding number of free parameters 
            // contained in each block, as listed in the input prior file
            fullPathOutputHyperParameters = outputPathPrefix + "hyperParameters.txt";
            ArrayXd priorIdentifiers = hyperParameters.col(4);
            ArrayXi NparametersPerPriorBlock(1);
            ArrayXi priorBlockIdentifier(1);
            int currentIdentifier = priorIdentifiers[0];
            int NpriorTypes = 1;                                    // At least one prior type must be included in the computation
            NparametersPerPriorBlock << 1;
            priorBlockIdentifier << currentIdentifier;

            for (int i = 1; i < priorIdentifiers.size(); i++)
            {
                if (priorIdentifiers[i] != currentIdentifier)
                {
                    currentIdentifier = priorIdentifiers[i];
                    NpriorTypes++;
                    NparametersPerPriorBlock.conservativeResize(NparametersPerPriorBlock.size() + 1);
                    NparametersPerPriorBlock[NparametersPerPriorBlock.size() - 1] = 1;
                    priorBlockIdentifier.conservativeResize(priorBlockIdentifier.size() + 1);
                    priorBlockIdentifier[priorBlockIdentifier.size() - 1] = currentIdentifier;
                }
                else
                {
                    NparametersPerPriorBlock[NparametersPerPriorBlock.size() - 1]++;
                }
            }

            cout << "------------------------------------------------------- " << endl;
            cout << " A total of " << NpriorTypes << " blocks of different" << endl; 
            cout << " prior types have been found." << endl;
            cout << "------------------------------------------------------- " << endl;

            vector<Prior*> ptrPriors(NpriorTypes);
            int segmentFirstIndex = 0;
    
            // Define prior hyper parameters for each block of a given prior type
            for (int i = 0; i < NpriorTypes; i++)
            {

                if (priorBlockIdentifier[i] == 0)
                {
                    // Uniform prior
                    ArrayXd hyperParametersMinima = hyperParameters.col(0).segment(segmentFirstIndex, NparametersPerPriorBlock[i]);
                    ArrayXd hyperParametersMaxima = hyperParameters.col(1).segment(segmentFirstIndex, NparametersPerPriorBlock[i]);

                    ptrPriors[i] = new UniformPrior(hyperParametersMinima, hyperParametersMaxima);
                }

                if (priorBlockIdentifier[i] == 1)
                {
                    // Gaussian prior
                    ArrayXd hyperParametersMean = hyperParameters.col(0).segment(segmentFirstIndex, NparametersPerPriorBlock[i]);
                    ArrayXd hyperParametersSDV = hyperParameters.col(1).segment(segmentFirstIndex, NparametersPerPriorBlock[i]);

                    ptrPriors[i] = new NormalPrior(hyperParametersMean, hyperParametersSDV);;
                }

                if (priorBlockIdentifier[i] == 2)
                {
                    // Super Gaussian prior
                    ArrayXd hyperParametersMean = hyperParameters.col(0).segment(segmentFirstIndex, NparametersPerPriorBlock[i]);
                    ArrayXd hyperParametersSDV = hyperParameters.col(1).segment(segmentFirstIndex, NparametersPerPriorBlock[i]);
                    ArrayXd hyperParametersWOP = hyperParameters.col(2).segment(segmentFirstIndex, NparametersPerPriorBlock[i]);

                    ptrPriors[i] = new SuperGaussianPrior(hyperParametersMean, hyperParametersSDV, hyperParametersWOP);  
                }

                if (priorBlockIdentifier[i] == 3)
                {
                    // Grid Uniform prior
                    ArrayXd hyperParametersStartingCoordinate = hyperParameters.col(0).segment(segmentFirstIndex, NparametersPerPriorBlock[i]);
                    ArrayXd hyperParametersEndingCoordinate = hyperParameters.col(1).segment(segmentFirstIndex, NparametersPerPriorBlock[i]);
                    ArrayXd hyperParametersNgridPoints = hyperParameters.col(2).segment(segmentFirstIndex, NparametersPerPriorBlock[i]);
                    ArrayXd hyperParametersTolerance = hyperParameters.col(3).segment(segmentFirstIndex, NparametersPerPriorBlock[i]);
                    
                    ptrPriors[i] = new GridUniformPrior(hyperParametersStartingCoordinate, hyperParametersEndingCoordinate, hyperParametersNgridPoints, hyperParametersTolerance);
                }
                
                segmentFirstIndex += NparametersPerPriorBlock[i];
            }

            if (writePriorHyperParametersToFile)
            {
                ofstream outputFile;
                File::openOutputFile(outputFile, fullPathOutputHyperParameters);
                outputFile << "# Hyper parameters used for setting up mixed prior distributions (Columns 1-4)." << endl;
                outputFile << "# Last column contains the prior type for each free parameter:" << endl;
                outputFile << "# 0 for Uniform, 1 for Gaussian, 2 for Super Gaussian, 3 for Grid Uniform." << endl;
                File::arrayXXdToFile(outputFile, hyperParameters);
                outputFile.close();
            }
           
            return ptrPriors;
        }
        else
        {
            cerr << " Wrong number of input prior hyper-parameters." << endl;
            cerr << " If only uniform priors are to be set, then a two-column file can be used." << endl;
            cerr << " If mixed prior types are to be used instead, then a 5-column file must be" << endl;
            cerr << " supplied, containing in the last column the number to identify the prior" << endl;
            cerr << " type of each free parameter (0 for Uniform, 1 for Gaussian, 2 for Super " << endl;
            cerr << " Gaussian, 3 for Grid Uniform)." << endl;
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        cout << "------------------------------------------------------- " << endl;
        cout << " Reading uniform prior file " + inputFileName << endl;
        cout << "------------------------------------------------------- " << endl;
        cout << endl;

        fullPathOutputHyperParameters = outputPathPrefix + "hyperParametersUniform.txt";
        ArrayXd hyperParametersMinima = hyperParameters.col(0);
        ArrayXd hyperParametersMaxima = hyperParameters.col(1);

        // Uniform Prior
        int NpriorTypes = 1;                                    // Only one prior type is included in the computation
        vector<Prior*> ptrPriors(NpriorTypes);
    
        ptrPriors[0] = new UniformPrior(hyperParametersMinima, hyperParametersMaxima);

        if (writePriorHyperParametersToFile)
        {
            ofstream outputFile;
            File::openOutputFile(outputFile, fullPathOutputHyperParameters);
            outputFile << "# Hyper parameters used for setting up uniform prior distributions (Col 1: minima, Col 2: maxima)." << endl;
            File::arrayXXdToFile(outputFile, hyperParameters);
            outputFile.close();
        }

        return ptrPriors;
    }
    
}
