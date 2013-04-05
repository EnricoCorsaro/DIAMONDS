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
: Prior(mean.size()),
  engine(time(0)),  
  mean(mean),
  standardDeviation(standardDeviation)
{
    assert (mean.size() == standardDeviation.size());
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
} // END NormalPrior::getNormalPrior()










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
    

} // END NormalPrior::drawWithConstraint()


