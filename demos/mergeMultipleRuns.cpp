//
// Compile with: clang++ -o merger mergeMultipleRuns.cpp -L../build/ -I ../include/ -l diamonds -stdlib=libc++ -std=c++11 -Wno-deprecated-register
// 

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Functions.h"
#include "File.h"
#include "EuclideanMetric.h"
#include "Prior.h"
#include "Results.h"
#include "ZeroModel.h"
#include "ZeroLikelihood.h"
#include "ZeroClusterer.h"
#include "ZeroSampler.h"
#include "ZeroPrior.h"

int main(int argc, char *argv[])
{
    // Check number of arguments for main function
    
    if (argc != 4)
    {
        cerr << "Usage: ./merger <KIC ID> <initial run number> <tot number of runs>" << endl;
        exit(EXIT_FAILURE);
    }

    // ---------------------------
    // ----- Read input data -----
    // ---------------------------

    string preName = "TO BE SET BY THE USER";
    string likelihoodFileName = preName + "logLikelihood.txt";
    string parameterFileName = preName + "parameter";
    string configurationFileName = preName + "configuringParameters.txt";

    unsigned long Npoints = 0;
    unsigned long NconfiguringParameters = 0;
    int Ncols = 0;
    string KIC_ID(argv[1]);
    string inputDirName = "TO BE SET BY THE USER";
    string outputDirName = inputDirName;          
    string initialProcessReference(argv[2]);
    int initialProcessNumber = atoi(initialProcessReference.c_str());
    string numberOfProcesses(argv[3]);      
    int Nprocesses = atoi(numberOfProcesses.c_str());   // The total number of processes to be merged (starting from 0)


    ArrayXd totalLogLikelihood;
    ArrayXXd totalParameterValues;
    int NfiguresProcesses = 2;

    cerr << "------------------------------------------------" << endl;
    cerr << " Total processes to be merged: " << Nprocesses << endl;
    cerr << " Starting process: " << initialProcessNumber << endl;
    cerr << "------------------------------------------------" << endl;
    cerr << endl;


    // Read dimensions of the inference analysis from one of the stored results

    ArrayXXd data;
    ostringstream numberString;
    numberString << setfill('0') << setw(NfiguresProcesses) << 0;
    string inputFileName3 = inputDirName + numberString.str() + "/" + configurationFileName;
    ifstream inputFile;
    File::openInputFile(inputFile, inputFileName3);
    File::sniffFile(inputFile, Npoints, Ncols);
    data = File::arrayXXdFromFile(inputFile, Npoints, Ncols);
    inputFile.close();

    int Ndimensions = data(0,0);
    int NfiguresDimensions = ceil(log10(Ndimensions*100));


    // Read input data from different processes and collect them
    
    int totalNpoints = 0;
    int totalNlive = 0;
    ArrayXXd data1;
    ArrayXXd data2;

    for (int i=initialProcessNumber; i < initialProcessNumber+Nprocesses; i++)
    {
        // Read log likelihood values
       
        ostringstream numberString1;
        numberString1 << setfill('0') << setw(NfiguresProcesses) << i;
        string inputFileName1 = inputDirName + numberString1.str() + "/" + likelihoodFileName;
        File::openInputFile(inputFile, inputFileName1);
        File::sniffFile(inputFile, Npoints, Ncols);
        data1 = File::arrayXXdFromFile(inputFile, Npoints, Ncols);
        inputFile.close();
        
        totalLogLikelihood.conservativeResize(totalNpoints + Npoints);
        totalLogLikelihood.segment(totalNpoints, Npoints) = data1;


        // Read configuring parameters to compute total number of live points

        string inputFileName3 = inputDirName + numberString1.str() + "/" + configurationFileName;
        File::openInputFile(inputFile, inputFileName3);
        File::sniffFile(inputFile, NconfiguringParameters, Ncols);
        data = File::arrayXXdFromFile(inputFile, NconfiguringParameters, Ncols);
        inputFile.close();

        int Nlive = data(1,0);
        totalNlive += Nlive;


        // Read parameter values for each dimension and collect them

        totalParameterValues.conservativeResize(Ndimensions, totalNpoints + Npoints);
        
        for (int j=0; j < Ndimensions; j++)
        {
            ostringstream numberString2;
            numberString2 << setfill('0') << setw(NfiguresDimensions) << j;
            string inputFileName2 = inputDirName + numberString1.str() + "/" + parameterFileName + numberString2.str() + ".txt";
            File::openInputFile(inputFile, inputFileName2);
            File::sniffFile(inputFile, Npoints, Ncols);
            data2 = File::arrayXXdFromFile(inputFile, Npoints, Ncols);
            inputFile.close();

            totalParameterValues.row(j).segment(totalNpoints, Npoints) = data2.transpose();
        }

        totalNpoints += Npoints;
    }
  

    // Sort total arrays by increasing likelihood and sort parameter values accordingly

    for (int i=0; i < Ndimensions; i++)
    {
        ArrayXd parameterValues(totalNpoints);
        parameterValues = totalParameterValues.row(i);
        Functions::topDownMergeSort(totalLogLikelihood, parameterValues);
        totalParameterValues.row(i) = parameterValues;
    }

    
    // Compute weights according to trapezoidal rule
    
    ArrayXd totalLogWeight(totalNpoints);
    double reductionFactor = exp(-1.0/totalNlive);

    totalLogWeight(0) = log(0.5) + log(2 - reductionFactor - reductionFactor*reductionFactor);         // First boundary condition
    totalLogWeight(totalNpoints-1) = log(0.5) + log(pow(reductionFactor,totalNpoints-2) - pow(reductionFactor,totalNpoints));      // Second boundary condition

    for (int i=1; i < totalNpoints - 1; i++)
    {
        totalLogWeight(i) = log(0.5) + log(pow(reductionFactor,i-1) - pow(reductionFactor,i+1));
    }
    

    // Compute final total evidence and its error bar

    double logEvidence = numeric_limits<double>::lowest();
    double informationGain = 0.0; 


    // Update the evidence and the information Gain
    
    for (int i=0; i < totalNpoints; ++i)
    {
        double logEvidenceContributionNew = totalLogWeight(i) + totalLogLikelihood(i);    
        double logEvidenceNew = Functions::logExpSum(logEvidence, logEvidenceContributionNew);
        informationGain = exp(logEvidenceContributionNew - logEvidenceNew) * totalLogLikelihood(i) 
                    + exp(logEvidence - logEvidenceNew) * (informationGain + logEvidence) 
                    - logEvidenceNew;
        logEvidence = logEvidenceNew;
    }


    // Compute Skilling's error on the log(Evidence)
    
    double logEvidenceError = sqrt(fabs(informationGain)/totalNlive); 



    // Set up NestedSampler with current information
    
    ArrayXd covariates;
    ArrayXd observations;

    ZeroModel model(covariates);
    ZeroLikelihood likelihood(observations, model);

    EuclideanMetric metric;
    ZeroClusterer clusterer(metric);
    
    bool printOnTheScreen = true;                   // Print results on the screen
    int initialNobjects = 500;                      // Maximum (and initial) number of live points evolving within the nested sampling process. 
    int minNobjects = 500;                          // Minimum number of live points allowed in the computation
    
    vector<Prior*> ptrPriors(1);
    ZeroPrior zeroPrior(1);
    ptrPriors[0] = &zeroPrior;

    ZeroSampler nestedSampler(printOnTheScreen, initialNobjects, minNobjects, ptrPriors, likelihood, metric, clusterer);
    
    nestedSampler.setLogEvidence(logEvidence);
    nestedSampler.setLogEvidenceError(logEvidenceError);
    nestedSampler.setInformationGain(informationGain);
    nestedSampler.setPosteriorSample(totalParameterValues);
    nestedSampler.setLogLikelihoodOfPosteriorSample(totalLogLikelihood);
    nestedSampler.setLogWeightOfPosteriorSample(totalLogWeight);
    nestedSampler.setOutputPathPrefix(outputDirName);
    

    string outputFileName = "mergedConfiguringParameters.txt";
    string fullPath = outputDirName + outputFileName;
    File::openOutputFile(nestedSampler.outputFile, fullPath);
    nestedSampler.outputFile << "# List of configuring parameters deriving from the merging of multiple runs." << endl;
    nestedSampler.outputFile << "# Row #1: Ndimensions" << endl;
    nestedSampler.outputFile << "# Row #2: Nprocesses" << endl;
    nestedSampler.outputFile << "# Row #3: Total Nobjects" << endl;
    nestedSampler.outputFile << "# Row #4: Total Niterations" << endl;
    nestedSampler.outputFile << Ndimensions << endl;
    nestedSampler.outputFile << Nprocesses << endl;
    nestedSampler.outputFile << totalNlive << endl;
    nestedSampler.outputFile << totalNpoints << endl;

    nestedSampler.outputFile.close();

    Results results(nestedSampler);
    results.writeParametersToFile("parameter");
    results.writeLogLikelihoodToFile("logLikelihood.txt");
    results.writeLogWeightsToFile("logWeight.txt");
    results.writeEvidenceInformationToFile("evidenceInformation.txt");
    results.writePosteriorProbabilityToFile("posteriorDistribution.txt");

    double credibleLevel = 68.3;
    bool writeMarginalDistributionToFile = true;
    results.writeParametersSummaryToFile("parameterSummary.txt", credibleLevel, writeMarginalDistributionToFile);

    cerr << "------------------------------------------------" << endl;
    cerr << " Merging complete." << endl;
    cerr << "------------------------------------------------" << endl;

    return EXIT_SUCCESS;
}
