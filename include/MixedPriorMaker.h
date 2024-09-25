// Namespace for preparation of mixed prior distributions.
// Enrico Corsaro @ OACT - 30 July 2024
// e-mail: emncorsaro@gmail.com
// Header file "MixedPriorMaker.h"
// Implementation contained in "MixedPriorMaker.cpp"

#ifndef MIXEDPRIORMAKER_H
#define MIXEDPRIORMAKER_H

#include "Functions.h"
#include "UniformPrior.h"
#include "NormalPrior.h"
#include "SuperGaussianPrior.h"
#include "GridUniformPrior.h"

namespace MixedPriorMaker
{
    vector<Prior*> prepareDistributions(string inputFileName, string outputPathPrefix, unsigned long &Ndimensions, 
                                             bool writePriorHyperParametersToFile);

}

#endif