#include "NormalPrior.h"



// NormalPrior:NormalPrior()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
//      mean: array containing mean values for setting 
//      centroids of the multi-dimensional normal prior. 
//      standardDeviation: array containing SDV values of the
//      multi-dimensional normal prior.
//

NormalPrior::NormalPrior(const RefArrayXd mean, const RefArrayXd standardDeviation)
: Prior(mean.size(),false),
  mean(mean),
  standardDeviation(standardDeviation)
{
    assert(mean.size() == standardDeviation.size());
    normalDistributionVector.resize(mean.size());

    for (int i = 0; i < mean.size(); i++)
    {
        normal_distribution<double> normal(mean(i),standardDeviation(i));
        normalDistributionVector[i] = normal;
    }

    normalizingFactor = (1./(standardDeviation*sqrt(2*Functions::PI))).prod();

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












// NormalPrior::getNormalizingFactor()
//
// PURPOSE: 
//      Get the private data member normalizingFactor.
//
// OUTPUT:
//      An integer containing the normalization factor of the normal prior.

double NormalPrior::getNormalizingFactor()
{
    return normalizingFactor;
}











// NormalPrior::draw()
//
// PURPOSE:
//      Draw a sample of parameters values from a normal prior
//      distribution. The parameters are in number Ndimensions
//      and contain Nobjects values each.
//
// INPUT:
//      nestedSampleOfParameters: two-dimensional Eigen Array to contain 
//      the resulting parameters values.
//      Nobjects: integer containing the number of objects to be drawn.
//
// OUTPUT:
//      void
//

void NormalPrior::draw(RefArrayXXd nestedSampleOfParameters, const int Nobjects)
{
    // Normal sampling over parameters intervals
    
    for (int i = 0; i < Ndimensions; i++)
    {
        for (int j = 0; j < Nobjects; j++)
        {
            nestedSampleOfParameters(i,j) = normalDistributionVector[i](engine);
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
//      nestedSampleOfParameters: one-dimensional array containing the set of 
//      parameters values to be updated.
//      likelihood: an object to compute the corresponding likelihood value.
//
// OUTPUT:
//      void
//
// NOTE:
//      nestedSampleOfParameters refers to the worst object identified in the nested
//      sampling loop. Thus, the array contains Ndimensions elements.
//

void NormalPrior::drawWithConstraint(RefArrayXd nestedSampleOfParameters, Likelihood &likelihood)
{
    double logLikelihood;
    double logLikelihoodConstraint = likelihood.logValue(nestedSampleOfParameters);


    // Normal sampling to find new parameter with logLikelihood > logLikelihoodConstraint
    
    do
    {
        for (int i = 0; i < Ndimensions; i++)
            {
                nestedSampleOfParameters(i) = normalDistributionVector[i](engine);
            }
    
        logLikelihood = likelihood.logValue(nestedSampleOfParameters);
    }
    while (logLikelihood <= logLikelihoodConstraint);
    

} 













// NormalPrior::pointIsRejected()
//
// PURPUSE:
//      Evaluates whether input point coordinates satisfy prior conditions.
//
// INPUT:
//      drawnSampleOfParameters: an Eigen Array of size (Ndimensions, 2)
//      containing a sample of coordinates for one object to be verified (column 0)
//      and for a reference object used in the sampling process (column 1).
//
// OUTPUT:
//      A bool variable declaring whether the point has to be rejected (true)
//      or accepted (false)
//

bool NormalPrior::pointIsRejected(RefArrayXXd drawnSampleOfParameters)
{
    assert(drawnSampleOfParameters.cols() == 2);
    assert(drawnSampleOfParameters.rows() == Ndimensions);
    
    bool pointIsRejected = false;

    
    // Compute density-related values for each point

    double weight1;
    double weight2;

    weight1 = ((drawnSampleOfParameters.col(0) - mean)/standardDeviation).square().sum();
    weight2 = ((drawnSampleOfParameters.col(1) - mean)/standardDeviation).square().sum();

    if (weight1 > weight2)          // If prior density of drawn point is lower than reference point
        pointIsRejected = true;     // then drawn point is not accepeted
                                                        
    return pointIsRejected;
}

