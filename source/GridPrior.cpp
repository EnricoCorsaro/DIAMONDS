#include "GridPrior.h"



// GridPrior::GridPrior()
//
// PURPOSE: 
//      Derived class constructor. Builds a grid prior with uniform sampling
//      only inside the steps of the grid.
//
// INPUT:
//      width:              an array containing the width of the grid steps. 
//      separation:         an array containing the separation of the grid steps. It has to be a number greater than the width itself.
//      startingCoordinate: an array containing the starting point of the grid prior in each coordinate.
//      Nsteps:             an array of integers containing the number of steps in each dimension.
//
GridPrior::GridPrior(const RefArrayXd width, const RefArrayXd separation, const RefArrayXd startingCoordinate,
                     const RefArrayXd Nsteps)
: Prior(width.size()),
  uniform(0.0,1.0),
  width(width),
  separation(separation),
  startingCoordinate(startingCoordinate),
  Nsteps(Nsteps)
{
    assert (width.size() == separation.size());
    assert (startingCoordinate.size() == width.size());

    if ( (width <= 0.0).any() || (separation <= width).any() )
    {
        cerr << "Grid Prior hyper parameters are not correctly typeset." << endl;
        exit(EXIT_FAILURE);
    }

    uniformIntegerVector.resize(width.size());

    for (int i = 0; i < width.size(); i++)
    {
        uniform_int_distribution<int> uniformInteger(0, Nsteps(i)-1);
        uniformIntegerVector[i] = uniformInteger;
    }
}









// GridPrior:~GridPrior()
//
// PURPOSE: 
//      Derived class destructor.
//

GridPrior::~GridPrior()
{

}









// GridPrior::getWidth()
//
// PURPOSE:
//      Get the private data member width.
//
// OUTPUT:
//      An array containing the width of the grid steps.

ArrayXd GridPrior::getWidth()
{
    return width;    
}











// GridPrior::getSeparation()
//
// PURPOSE:
//      Get the private data member separation.
//
// OUTPUT:
//      An array containing the separation of the grid steps. It has to be a number greater than the width itself.
//

ArrayXd GridPrior::getSeparation()
{
    return separation;    
}










// GridPrior::getStartingCoordinate()
//
// PURPOSE:
//      Get the private data member startingCoordinate.
//
// OUTPUT:
//      An array containing the starting point of the grid prior in each coordinate.
//

ArrayXd GridPrior::getStartingCoordinate()
{
    return startingCoordinate;    
}











// GridPrior::getNsteps()
//
// PURPOSE:
//      Get the private data member Nsteps.
//
// OUTPUT:
//      An array of integers containing the number of steps in each dimension.
//

ArrayXd GridPrior::getNsteps()
{
    return Nsteps;    
}










// GridPrior::logDensity()
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

double GridPrior::logDensity(RefArrayXd const x, const bool includeConstantTerm)
{
    double logDens;
    ArrayXd modulo(width.size());

    for (int i=0; i < width.size(); ++i)
    {
        modulo(i) = fmod(x(i) - startingCoordinate(i), separation(i));
    }

    if ((modulo > width).any() || (x < startingCoordinate).any() || (x > startingCoordinate + Nsteps*separation).any())
    {
        // The point x falls out of the grid steps. In this case
        // the density is zero, and thus the log(density) is -infinity. The value
        // of the latter is defined in the parent class.

        logDens = minusInfinity;
        return logDens;
    }
    else
    {
        // The point falls inside the distribution's boundaries.

        logDens = 0.0;
    }

    if (includeConstantTerm)
    {
        logDens += (-1.0) * (Nsteps*width).log().sum(); 
    }

    return logDens;
}













// GridPrior::drawnPointIsAccepted()
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

bool GridPrior::drawnPointIsAccepted(RefArrayXd const drawnPoint)
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












// GridPrior::draw()
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

void GridPrior::draw(RefArrayXXd drawnSample)
{
    int Npoints = drawnSample.cols();


    // Grid sampling over all free parameters and points

    for (int i = 0; i < Ndimensions; i++)
    {
        for (int j = 0; j < Npoints; j++)
        {
            drawnSample(i,j) = width(i)*uniform(engine) + uniformIntegerVector[i](engine)*separation(i) + startingCoordinate(i);
        }
    }

}







// GridPrior::drawWithConstraint()
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

void GridPrior::drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood)
{
    double logLikelihood;
    double logLikelihoodConstraint = likelihood.logValue(drawnPoint);


    // Grid sampling to find new parameter with logLikelihood > logLikelihoodConstraint
    
    do
    {
        for (int i = 0; i < Ndimensions; i++)
            {
                drawnPoint(i) = width(i)*uniform(engine) + uniformIntegerVector[i](engine)*separation(i) + startingCoordinate(i);
            }
    
        logLikelihood = likelihood.logValue(drawnPoint);
    }
    while (logLikelihood <= logLikelihoodConstraint);
    

}












// GridPrior::writeHyperParametersToFile()
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

void GridPrior::writeHyperParametersToFile(string fullPath)
{
    ArrayXXd hyperParameters(width.size(),4);
    hyperParameters.col(0) = width;
    hyperParameters.col(1) = separation;
    hyperParameters.col(2) = startingCoordinate;
    hyperParameters.col(3) = Nsteps*1.0;

    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
    outputFile << "# Hyper parameters used for setting up grid priors." << endl;
    outputFile << "# Each line corresponds to a different free parameter (coordinate)." << endl;
    outputFile << "# Column #1: Width (of the steps)" << endl;
    outputFile << "# Column #2: Separation (between different steps)" << endl;
    outputFile << "# Column #3: Starting coordinate" << endl;
    outputFile << "# Column #4: Number of steps in the grid" << endl;
    File::arrayXXdToFile(outputFile, hyperParameters);
    outputFile.close();
}
