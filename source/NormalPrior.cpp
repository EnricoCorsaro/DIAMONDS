#include "NormalPrior.h"



// NormalPrior:NormalPrior()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
//      mean: array containing mean values for setting 
//            centroids of the multi-dimensional normal prior. 
//      standardDeviation: array containing SDV values of the
//                         multi-dimensional normal prior.
//

NormalPrior::NormalPrior(const RefArrayXd mean, const RefArrayXd standardDeviation)
: Prior(mean.size()),
  mean(mean),
  standardDeviation(standardDeviation)
{
    assert(mean.size() == standardDeviation.size());
    assert(standardDeviation.any() >= 0);
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

double NormalPrior::logDensity(RefArrayXd x, const bool includeConstantTerm)
{
    double logDens = -0.5 * ((x - mean)/standardDeviation).square().sum();

    if (includeConstantTerm)
    {
        logDens += -Ndimensions/2.0 * log(2.*Functions::PI) - 0.5 * standardDeviation.log().sum();
    }

    return logDens;
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








