#include "NormalPrior.h"



// NormalPrior:NormalPrior()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
//      mean:               array containing mean values for setting 
//                          centroids of the multi-dimensional normal prior. 
//      standardDeviation:  array containing SDV values of the
//                          multi-dimensional normal prior.
//

NormalPrior::NormalPrior(RefArrayXd const mean, RefArrayXd const standardDeviation)
: Prior(mean.size()),
  uniform(0.0,1.0),
  mean(mean),
  standardDeviation(standardDeviation)
{
    assert(mean.size() == standardDeviation.size());

    if ( (standardDeviation <= 0.0).any() )
    {
        cerr << "Normal Prior hyper parameters are not correctly typeset." << endl;
        exit(EXIT_FAILURE);
    }

    normalDistributionVector.resize(mean.size());

    for (int i = 0; i < mean.size(); i++)
    {
        normal_distribution<double> normal(mean(i),standardDeviation(i));
        normalDistributionVector[i] = normal;
    }
}











// NormalPrior:~NormalPrior()
//
// PURPOSE: 
//      Derived class destructor.
//

NormalPrior::~NormalPrior()
{

}












// NormalPrior::getMean()
//
// PURPOSE:
//      Get the private data member mean.
//
// OUTPUT:
//      An array containing the mean values of the normal prior
//      for each of the free parameters.
//

ArrayXd NormalPrior::getMean()
{
    return mean;    
}











// NormalPrior::getStandardDeviation()
//
// PURPOSE:
//      Get the private data member standardDeviation.
//
// OUTPUT:
//      An array containing the SDV values of the normal prior
//      for each of the free parameters.
//

ArrayXd NormalPrior::getStandardDeviation()
{
    return standardDeviation;    
}










// NormalPrior::logDensity()
//
// PURPOSE:
//      Compute the logarithm of the probability density distribution evaluated in 'x'.
//
// INPUT: 
//      x:                   Point in which the log(pdf) should be evaluated.
//      includeConstantTerm: If true : compute the exact log(density), 
//                           If false: ignore the constant terms (with factors of pi, 2, etc.)
//
// OUTPUT:
//      Natural logarithm of the probability density evaluation in x.
//

double NormalPrior::logDensity(RefArrayXd const x, const bool includeConstantTerm)
{
    assert(x.size() == Ndimensions);

    double logDens = -0.5 * ((x - mean)/standardDeviation).square().sum();

    if (includeConstantTerm)
    {
        logDens += -Ndimensions/2.0 * log(2.*Functions::PI) - standardDeviation.log().sum();
    }

    return logDens;
}













// NormalPrior::drawnPointIsAccepted()
//
// PURPOSE:
//      Checks wheter a new drawn point to be verified is accepted 
//      according to the prior distribution.
//
// INPUT: 
//      drawnPoint:     an Eigen array containing the coordinates 
//                      of the new drawn point to be verified
//
// OUTPUT:
//      Returns true if the new point is accepted, false if not.
//

bool NormalPrior::drawnPointIsAccepted(RefArrayXd const drawnPoint)
{
    // Evaluate the distribution density value in a normalized scale 

    double normalizedDensity = exp(logDensity(drawnPoint));

    
    // Compute a reference density from a uniform distribution between 0 and 1
    
    double referenceNormalizedDensity = uniform(engine);
    

    // Compare the two densities and accept the point only if the drawn density is larger than the reference one

    bool pointIsAccepted = normalizedDensity > referenceNormalizedDensity;

    return pointIsAccepted;
}












// NormalPrior::draw()
//
// PURPOSE:
//      Draw a sample of parameters values from a normal prior
//      distribution. The parameters are in number Ndimensions
//      and contain Npoints values each.
//
// INPUT:
//      drawnSample:    two-dimensional Eigen Array to contain 
//                      the resulting parameters values.
//
// OUTPUT:
//      void
//

void NormalPrior::draw(RefArrayXXd drawnSample)
{ 
    int Npoints = drawnSample.cols();
 
    
    // Normal sampling over all free parameters and points 

    for (int i = 0; i < Ndimensions; i++)
    {
        for (int j = 0; j < Npoints; j++)
        {
            drawnSample(i,j) = normalDistributionVector[i](engine);
        }
    }

}







// NormalPrior::drawWithConstraint()
//
// PURPOSE: 
//      Replace an old set of parameters values with a new one
//      having higher likelihood value.
//
// INPUT:
//      drawnPoint: one-dimensional array containing the set of parameters
//                  values to be updated. Initially it should contain the coordinates
//                  of the point with the worst likelihood.
//      likelihood: an object to compute the likelihood values.
//
// OUTPUT:
//      void

void NormalPrior::drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood)
{
    double logLikelihood;
    double logLikelihoodConstraint = likelihood.logValue(drawnPoint);


    // Normal sampling to find new parameter with logLikelihood > logLikelihoodConstraint
    
    do
    {
        for (int i = 0; i < Ndimensions; i++)
        {
            drawnPoint(i) = normalDistributionVector[i](engine);
        }
    
        logLikelihood = likelihood.logValue(drawnPoint);
    }
    while (logLikelihood <= logLikelihoodConstraint);
} 









// NormalPrior::writeHyperParametersToFile()
//
// PURPOSE: 
//      Store prior hyper parameters in an output ASCII file.
//
// INPUT:
//      fullPath:     a string containing the full filename for the output ASCII file
//
// OUTPUT:
//      void
//

void NormalPrior::writeHyperParametersToFile(string fullPath)
{
    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
    outputFile << "# Hyper parameters used for setting up normal priors." << endl;
    outputFile << "# Each line corresponds to a different free parameter (coordinate)." << endl;
    outputFile << "# Column #1: Mean" << endl;
    outputFile << "# Column #2: Standard Deviation" << endl;
    File::twoArrayXdToFile(outputFile, mean, standardDeviation);
    outputFile.close();
}
