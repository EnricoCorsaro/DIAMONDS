#include "GridUniformPrior.h"



// GridUniformPrior::GridUniformPrior()
//
// PURPOSE: 
//      Derived class constructor. Builds a grid prior with uniform sampling
//      only inside the gridPoints of the grid.
//
// INPUT:
//      startingCoordinate: an array containing the starting point of the grid prior in each coordinate.
//      NgridPoints:        an array containing the number of grid points in each dimension (considered as an integer number).
//      separation:         an array containing the separation of the grid points. It has to be a number greater than the tolerance itself.
//      tolerance:          an array containing the tolerance of the grid points. This tolerance is an interval on each side of the grid point. 
//
GridUniformPrior::GridUniformPrior(const RefArrayXd startingCoordinate, const RefArrayXd NgridPoints, const RefArrayXd separation, 
                                   const RefArrayXd tolerance)
: Prior(tolerance.size()),
  uniform(-1.0,1.0),
  startingCoordinate(startingCoordinate),
  NgridPoints(NgridPoints),
  separation(separation),
  tolerance(tolerance)
{
    assert (tolerance.size() == separation.size());
    assert (startingCoordinate.size() == tolerance.size());

    if ( (tolerance <= 0.0).any() || (tolerance >= separation/2.).any() )
    {
        cerr << "Grid Prior hyper parameters are not correctly typeset." << endl;
        exit(EXIT_FAILURE);
    }

    uniformIntegerVector.resize(Ndimensions);

    for (int i = 0; i < Ndimensions; i++)
    {
        uniform_int_distribution<int> uniformInteger(0, NgridPoints(i)-1);
        uniformIntegerVector[i] = uniformInteger;
    }
}









// GridUniformPrior:~GridUniformPrior()
//
// PURPOSE: 
//      Derived class destructor.
//

GridUniformPrior::~GridUniformPrior()
{

}













// GridUniformPrior::getStartingCoordinate()
//
// PURPOSE:
//      Get the private data member startingCoordinate.
//
// OUTPUT:
//      An array containing the starting point of the grid prior in each coordinate.
//

ArrayXd GridUniformPrior::getStartingCoordinate()
{
    return startingCoordinate;    
}











// GridUniformPrior::getNgridPoints()
//
// PURPOSE:
//      Get the private data member NgridPoints.
//
// OUTPUT:
//      An array of integers containing the number of gridPoints in each dimension.
//

ArrayXd GridUniformPrior::getNgridPoints()
{
    return NgridPoints;    
}












// GridUniformPrior::getSeparation()
//
// PURPOSE:
//      Get the private data member separation.
//
// OUTPUT:
//      An array containing the separation of the grid gridPoints. It has to be a number greater than the tolerance itself.
//

ArrayXd GridUniformPrior::getSeparation()
{
    return separation;    
}












// GridUniformPrior::getTolerance()
//
// PURPOSE:
//      Get the private data member tolerance.
//
// OUTPUT:
//      An array containing the tolerance of the grid gridPoints.

ArrayXd GridUniformPrior::getTolerance()
{
    return tolerance;    
}












// GridUniformPrior::logDensity()
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

double GridUniformPrior::logDensity(RefArrayXd const x, const bool includeConstantTerm)
{
    double logDens;
    
    if ((x < (startingCoordinate - tolerance)).any() || (x > (startingCoordinate + (NgridPoints-1)*separation + tolerance)).any())
    {
        // The point x falls out of the grid of points (including tolerance on each side). In this case
        // the density is zero, and thus the log(density) is -infinity. The value
        // of the latter is defined in the parent class.

        logDens = minusInfinity;
        return logDens;
    }
    

    // Evaluate the distance of the given drawn point from the closest grid point 
    // to understand if its coordinates match the grid prior boundaries

    ArrayXd absoluteDistanceFromGridPoint(Ndimensions);

    for (int i = 0; i < Ndimensions; ++i)
    {
        absoluteDistanceFromGridPoint(i) = fmod(x(i) - startingCoordinate(i), separation(i));
        
        if (absoluteDistanceFromGridPoint(i) >= separation(i)/2.0)
            absoluteDistanceFromGridPoint(i) = separation(i) - absoluteDistanceFromGridPoint(i);

        if (absoluteDistanceFromGridPoint(i) >= tolerance(i))
        {
            logDens = minusInfinity;
            return logDens;
        }
    }
            
    // If you arrived until here, then the point falls within the boundaries of the distribution.

    logDens = 0.0;

    
    // Add the normalization factor of the integral if requested

    if (includeConstantTerm)
    {
        logDens += (-1.0) * (NgridPoints*tolerance*2).log().sum(); 
    }

    return logDens;
}













// GridUniformPrior::drawnPointIsAccepted()
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

bool GridUniformPrior::drawnPointIsAccepted(RefArrayXd const drawnPoint)
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












// GridUniformPrior::draw()
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

void GridUniformPrior::draw(RefArrayXXd drawnSample)
{
    int Npoints = drawnSample.cols();


    // Grid sampling over all free parameters and points

    for (int i = 0; i < Ndimensions; i++)
    {
        for (int j = 0; j < Npoints; j++)
        {
            drawnSample(i,j) = tolerance(i)*uniform(engine) + uniformIntegerVector[i](engine)*separation(i) + startingCoordinate(i);
        }
    }

}









// GridUniformPrior::drawWithConstraint()
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

void GridUniformPrior::drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood)
{
    double logLikelihood;
    double logLikelihoodConstraint = likelihood.logValue(drawnPoint);


    // Grid sampling to find new parameter with logLikelihood > logLikelihoodConstraint
    
    do
    {
        for (int i = 0; i < Ndimensions; i++)
            {
                drawnPoint(i) = tolerance(i)*uniform(engine) + uniformIntegerVector[i](engine)*separation(i) + startingCoordinate(i);
            }
    
        logLikelihood = likelihood.logValue(drawnPoint);
    }
    while (logLikelihood <= logLikelihoodConstraint);
    

}












// GridUniformPrior::writeHyperParametersToFile()
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

void GridUniformPrior::writeHyperParametersToFile(string fullPath)
{
    ArrayXXd hyperParameters(Ndimensions, 4);
    hyperParameters.col(0) = startingCoordinate;
    hyperParameters.col(1) = NgridPoints;
    hyperParameters.col(2) = tolerance;
    hyperParameters.col(3) = separation;

    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
    outputFile << "# Hyper parameters used for setting up grid priors." << endl;
    outputFile << "# Each line corresponds to a different free parameter (coordinate)." << endl;
    outputFile << "# Column #1: Tolerance (of the position of the grid points)" << endl;
    outputFile << "# Column #2: Separation (between different grid points)" << endl;
    outputFile << "# Column #3: Starting coordinate" << endl;
    outputFile << "# Column #4: Number of grid points in the grid" << endl;
    File::arrayXXdToFile(outputFile, hyperParameters);
    outputFile.close();
}
