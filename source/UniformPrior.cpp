#include "UniformPrior.h"



// UniformPrior::UniformPrior()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
//      minima:     array containing minimum values for setting 
//                  lower bounds of the free parameters. 
//      maxima:     array containing maximum values for setting
//                  upper bounds of the free parameters.
//

UniformPrior::UniformPrior(const RefArrayXd minima, const RefArrayXd maxima)
: Prior(minima.size()),
  uniform(0.0,1.0),
  minima(minima),
  maxima(maxima)
{
    assert (minima.size() == maxima.size());
    
    if ( (minima >= maxima).any() )
    {
        cerr << "Uniform Prior hyper parameters are not correctly typeset." << endl;
        exit(EXIT_FAILURE);
    }
}









// UniformPrior:~UniformPrior()
//
// PURPOSE: 
//      Derived class destructor.
//

UniformPrior::~UniformPrior()
{

}









// UniformPrior::getMinima()
//
// PURPOSE:
//      Get the private data member minima.
//
// OUTPUT:
//      An array containing the minimum values of the free parameters.
//

ArrayXd UniformPrior::getMinima()
{
    return minima;    
}











// UniformPrior::getMaxima()
//
// PURPOSE:
//      Get the private data member maxima.
//
// OUTPUT:
//      An array containing the maximum values of the free parameters.
//

ArrayXd UniformPrior::getMaxima()
{
    return maxima;    
}










// UniformPrior::logDensity()
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

double UniformPrior::logDensity(RefArrayXd const x, const bool includeConstantTerm)
{
    double logDens;

    if ((x < minima).any() | (x > maxima).any())
    {
        // The point x falls out of the distribution's boundaries. In this case
        // the density is zero, and thus the log(density) is -infinity. The value
        // of the latter is defined in the parent class.

        logDens = minusInfinity;
        return logDens;
    }
    else
    {
        // The point falls inside the boundaries of the distribution.

        logDens = 0.0;
    }

    if (includeConstantTerm)
    {
        logDens += (-1.0) * (maxima - minima).log().sum(); 
    }

    return logDens;
}













// UniformPrior::drawnPointIsAccepted()
//
// PURPOSE:
//      Checks whether a new drawn point to be verified is accepted 
//      according to the prior distribution.
//
// INPUT: 
//      drawnPoint:     an Eigen array containing the coordinates 
//                      of the new drawn point to be verified
//
// OUTPUT:
//      Returns true if the new point is accepted, false if not.
//

bool UniformPrior::drawnPointIsAccepted(RefArrayXd const drawnPoint)
{
    if (logDensity(drawnPoint) != minusInfinity)
    {
        return true;
    }
    else
    {
        return false;
    }
}












// UniformPrior::draw()
//
// PURPOSE:
//      Draw a sample of parameters values from a uniform prior
//      distributions. The parameters are in number Ndimensions
//      and contain Npoints values each.
//
// INPUT:
//      drawnSample:    two-dimensional Eigen Array to contain 
//                      the resulting parameters values.
//
// OUTPUT:
//      void
//

void UniformPrior::draw(RefArrayXXd drawnSample)
{
    int Npoints = drawnSample.cols();
 

    // Uniform sampling over all free parameters and points

    for (int i = 0; i < Ndimensions; i++)
    {
        for (int j = 0; j < Npoints; j++)
        {
            drawnSample(i,j) = uniform(engine)*(maxima(i)-minima(i)) + minima(i);
        }
    }

}







// UniformPrior::drawWithConstraint()
//
// PURPOSE: 
//      Replace an old set of parameters values with a new one
//      having higher likelihood value.
//
// INPUT:
//      drawnPoint:     one-dimensional Eigen Array containing the set of 
//                      parameters values to be updated.
//      likelihood:     an object to compute the corresponding likelihood value.
//
// OUTPUT:
//      void
//
// NOTE:
//      parameters refers to the worst object identified in the nested
//      sampling loop. Thus, the array contains Ndimensions elements.
//

void UniformPrior::drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood)
{
    double logLikelihood;
    double logLikelihoodConstraint = likelihood.logValue(drawnPoint);


    // Uniform sampling to find new parameter with logLikelihood > logLikelihoodConstraint
    
    do
    {
        for (int i = 0; i < Ndimensions; i++)
            {
                drawnPoint(i) = uniform(engine)*(maxima(i) - minima(i)) + minima(i);
            }
    
        logLikelihood = likelihood.logValue(drawnPoint);
    }
    while (logLikelihood <= logLikelihoodConstraint);
    

}












// UniformPrior::writeHyperParametersToFile()
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

void UniformPrior::writeHyperParametersToFile(string fullPath)
{
    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
    outputFile << "# Hyper parameters used for setting up uniform priors." << endl;
    outputFile << "# Each line corresponds to a different free parameter (coordinate)." << endl;
    outputFile << "# Column #1: Minima (lower boundaries)" << endl;
    outputFile << "# Column #2: Maxima (upper boundaries)" << endl;
    File::twoArrayXdToFile(outputFile, minima, maxima);
    outputFile.close();
}
