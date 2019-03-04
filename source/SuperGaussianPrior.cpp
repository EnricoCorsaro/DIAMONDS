#include "SuperGaussianPrior.h"



// SuperGaussianPrior:SuperGaussianPrior()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
//      center:             array containing values of the center for constructing 
//                          the multi-dimensional Super-Gaussian prior. 
//      sigma:              array containing the sigma values corresponding to the
//                          width of wings of the multi-dimensional Super-Gaussian prior.
//      widthOfPlateau:     array containing the widths of the plateau, i.e. the flat
//                          region of the multi-dimensional Super-Gaussian prior.
//

SuperGaussianPrior::SuperGaussianPrior(const RefArrayXd center, const RefArrayXd sigma, const RefArrayXd widthOfPlateau)
: Prior(center.size()),
  uniform(0.0,1.0),
  center(center),
  sigma(sigma),
  widthOfPlateau(widthOfPlateau),
  halfWidthOfPlateau(widthOfPlateau/2.0)
{
    assert(center.size() == sigma.size());
    assert(center.size() == widthOfPlateau.size());
    normalDistributionVector.resize(Ndimensions);

    if ( (sigma <= 0.0).any() || (widthOfPlateau <= 0.0).any() )
    {
        cerr << "Super Gaussian Prior hyper parameters are not correctly typeset." << endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < Ndimensions; i++)
    {
        normal_distribution<double> normal(center(i),sigma(i));
        normalDistributionVector[i] = normal;
    }


    // Compute the exact integral of the prior for the regions containing the Gaussian tails
    // Normalize areas of different prior regions (i.e. either tails or plateau) to total prior area

    tailsArea = sqrt(2*Functions::PI) * sigma;      
    totalArea = tailsArea + widthOfPlateau;
    plateauArea = widthOfPlateau/totalArea;
    tailsArea /= totalArea;
        
    normalizedAreas.resize(2,Ndimensions);
    normalizedAreas.row(0) = tailsArea;
    normalizedAreas.row(1) = plateauArea;
}











// SuperGaussianPrior:~SuperGaussianPrior()
//
// PURPOSE: 
//      Derived class destructor.
//

SuperGaussianPrior::~SuperGaussianPrior()
{

}












// SuperGaussianPrior::getCenter()
//
// PURPOSE:
//      Get the private data member center.
//
// OUTPUT:
//      An array containing the center values of the Super-Gaussian prior
//      for each of the free parameters.
//

ArrayXd SuperGaussianPrior::getCenter()
{
    return center;    
}











// SuperGaussianPrior::getSigma()
//
// PURPOSE:
//      Get the private data member sigma.
//
// OUTPUT:
//      An array containing the sigma values of the Super-Gaussian prior
//      for each of the free parameters.
//

ArrayXd SuperGaussianPrior::getSigma()
{
    return sigma;    
}










// SuperGaussianPrior::getWidthOfPlateau()
//
// PURPOSE:
//      Get the private data member widthOfPlateau.
//
// OUTPUT:
//      An array containing the width of the plateau of the Super-Gaussian prior
//      for each of the free parameters.
//

ArrayXd SuperGaussianPrior::getWidthOfPlateau()
{
    return widthOfPlateau;    
}













// SuperGaussianPrior::logDensity()
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

double SuperGaussianPrior::logDensity(RefArrayXd const x, const bool includeConstantTerm)
{

    // Determine in which region of the prior the point is falling
    
    ArrayXd position = x - center;
    double logDens = 0.0;
    
    for (int i = 0; i < Ndimensions; ++i)
    {
        
        if (fabs(position(i)) < halfWidthOfPlateau(i))
        {
            // If point is falling in the plateau region, give a constant prior having unitary amplitude
            
            continue;       // The total logDensity is uneffected by this case
        }
        else
        {
            // If point is falling in the Gaussian tails then compute the logDensity according to the non-normalized Gaussian function
            // having unitary maximum amplitude

            logDens += -0.5 * pow((fabs(position(i)) - halfWidthOfPlateau(i))/sigma(i),2);     // The total logDensity is updated
        }
    }


    if (includeConstantTerm)
    {
        logDens += -1.0 * (widthOfPlateau + sqrt(2.*Functions::PI)*sigma).log().sum();
    }

    return logDens;
}













// SuperGaussianPrior::drawnPointIsAccepted()
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

bool SuperGaussianPrior::drawnPointIsAccepted(RefArrayXd const drawnPoint)
{
    // Evaluate the distribution density value in a normalized scale 

    double normalizedDensity = exp(logDensity(drawnPoint));
   

    // Compute a reference density from a uniform distribution between 0 and 1
    
    double referenceNormalizedDensity = uniform(engine);
    

    // Compare the two densities and accept the point only if the drawn density is larger than the reference one

    bool pointIsAccepted = normalizedDensity > referenceNormalizedDensity;

    return pointIsAccepted;
}












// SuperGaussianPrior::draw()
//
// PURPOSE:
//      Draw a sample of parameters values from a Super-Gaussian prior
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

void SuperGaussianPrior::draw(RefArrayXXd drawnSample)
{
    int Npoints = drawnSample.cols();


    // Create random distribution of [0,1] possible integers to select starting region randomly

    uniform_int_distribution<int> uniform_int(0,1);

    
    // Sampling over all free parameters and points 

    for (int i = 0; i < Ndimensions; i++)
    {
        for (int j = 0; j < Npoints; j++)
        {
            // Select either tails or plateau of the distribution according to their prior mass

            int indexOfSelectedRegion = uniform_int(engine);
            double selectedNormalizedArea = normalizedAreas(indexOfSelectedRegion,i);
            double uniformNumber = uniform(engine);

            if (selectedNormalizedArea < uniformNumber)
            {
                selectedNormalizedArea = normalizedAreas(1 - indexOfSelectedRegion,i);      // If region is not selected, take the other one
            }

            double centralCoordinate = center(i);

            if (indexOfSelectedRegion == 0)
            {
                // In the case tails are selected sample the parameter according to the 
                // normal distribution centered at "center" and having SDV = "sigma"

                double normalDistributedCoordinate = normalDistributionVector[i](engine);
                
                while (normalDistributedCoordinate == centralCoordinate)
                {
                    // No sampling at values corresponding to the maximum of the normal distribution
                    // because these values are already accounted for in the plateau
                
                    if (normalDistributedCoordinate > center(i))
                    {
                        // In the case drawn point falls in the right-hand part, shift it to the right 
                        // by half of the width of the plateau

                        normalDistributedCoordinate += halfWidthOfPlateau(i);
                    }
                    else
                        if (normalDistributedCoordinate < center(i))
                        {
                            // In the case drawn point falls in the left-hand part, shift it to the left 
                            // by half of the width of the plateau
                    
                            normalDistributedCoordinate -= halfWidthOfPlateau(i);
                        }

                    drawnSample(i,j) = normalDistributedCoordinate;
                }
            }
            else
            {
                // In the case plateau are selected sample the parameter according to the 
                // uniform distribution centered at "center" and having width = "widthOfPlateau"

                double uniformDistributedCoordinate = uniform(engine)*widthOfPlateau(i) - halfWidthOfPlateau(i) + centralCoordinate; 
                
                drawnSample(i,j) = uniformDistributedCoordinate;
            }
        }
    }

}







// SuperGaussianPrior::drawWithConstraint()
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

void SuperGaussianPrior::drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood)
{
    double logLikelihood;
    double logLikelihoodConstraint = likelihood.logValue(drawnPoint);


    // Sampling to find new parameter with logLikelihood > logLikelihoodConstraint
    
    do
    {
        // Create random distribution of [0,1] possible integers to select starting region randomly

        uniform_int_distribution<int> uniform_int(0,1);

    
        // Sampling over all free parameters

        for (int i = 0; i < Ndimensions; i++)
        {
            // Select either tails or plateau of the distribution according to their prior mass

            int indexOfSelectedRegion = uniform_int(engine);
            double selectedNormalizedArea = normalizedAreas(indexOfSelectedRegion,i);
            double uniformNumber = uniform(engine);

            if (selectedNormalizedArea < uniformNumber)
            {
                selectedNormalizedArea = normalizedAreas(1 - indexOfSelectedRegion,i);      // If region is not selected, take the other one
            }

            double centralCoordinate = center(i);

            if (indexOfSelectedRegion == 0)
            {
                // In the case tails are selected sample the parameter according to the 
                // normal distribution centered at "center" and having SDV = "sigma"

                double normalDistributedCoordinate = normalDistributionVector[i](engine);
                
                while (normalDistributedCoordinate == centralCoordinate)
                {
                    // No sampling at values corresponding to the maximum of the normal distribution
                    // because these values are already accounted for in the plateau
                
                    if (normalDistributedCoordinate > center(i))
                    {
                        // In the case drawn point falls in the right-hand part, shift it to the right 
                        // by half of the width of the plateau

                        normalDistributedCoordinate += 0.5*widthOfPlateau(i);
                    }
                    else
                        if (normalDistributedCoordinate < center(i))
                        {
                            // In the case drawn point falls in the left-hand part, shift it to the left 
                            // by half of the width of the plateau
                    
                            normalDistributedCoordinate -= 0.5*widthOfPlateau(i);
                        }

                    drawnPoint(i) = normalDistributedCoordinate;
                }
            }
            else
            {
                // In the case plateau are selected sample the parameter according to the 
                // uniform distribution centered at "center" and having width = "widthOfPlateau"

                double uniformDistributedCoordinate = uniform(engine)*widthOfPlateau(i) - 0.5*widthOfPlateau(i) + centralCoordinate; 
                
                drawnPoint(i) = uniformDistributedCoordinate;
            }
        }


        logLikelihood = likelihood.logValue(drawnPoint);
    }
    while (logLikelihood <= logLikelihoodConstraint);
} 










// SuperGaussianPrior::writeHyperParametersToFile()
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

void SuperGaussianPrior::writeHyperParametersToFile(string fullPath)
{
    ArrayXXd hyperParameters(Ndimensions,3);
    hyperParameters.col(0) = center;
    hyperParameters.col(1) = sigma;
    hyperParameters.col(2) = widthOfPlateau;

    ofstream outputFile;
    File::openOutputFile(outputFile, fullPath);
    outputFile << "# Hyper parameters used for setting up super Gaussian priors." << endl;
    outputFile << "# Each line corresponds to a different free parameter (coordinate)." << endl;
    outputFile << "# Column #1: Center (center of the plateau region)" << endl;
    outputFile << "# Column #2: Sigma (SDV of the symmetric Gaussian tails)" << endl;
    outputFile << "# Column #3: Width of the plateau (total width, centered in \"center\")" << endl;
    File::arrayXXdToFile(outputFile, hyperParameters);
    outputFile.close();
}
