#include "UniformPrior.h"



// UniformPrior:UniformPrior()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
//      minima: array containing minimum values for setting 
//      lower bounds of the free parameters. 
//      maxima: array containing maximum values for setting
//      upper bounds of the free parameters.
//

UniformPrior::UniformPrior(const RefArrayXd minima, const RefArrayXd maxima)
: Prior(minima.size(), true),
  uniform(0.0,1.0),
  minima(minima),
  maxima(maxima)
{
    assert (minima.size() == maxima.size());
    normalizingFactor = (1./(maxima - minima)).prod();
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












// UniformPrior::getNormalizingFactor()
//
// PURPOSE: 
//      Get the private data member uniformFactor.
//
// OUTPUT:
//      An integer containing the normalization factor of the uniform prior.

double UniformPrior::getNormalizingFactor()
{
    return normalizingFactor;
}












// UniformPrior::draw()
//
// PURPOSE:
//      Draw a sample of parameters values from a uniform prior
//      distributions. The parameters are in number Ndimensions
//      and contain Nobjects values each.
//
// INPUT:
//      nestedSampleOfParameters: two-dimensional Eigen Array to contain 
//      the resulting parameters values.
//
// OUTPUT:
//      void
//

void UniformPrior::draw(RefArrayXXd nestedSampleOfParameters)
{
    // Uniform sampling over parameters intervals
    
    for (int i = 0; i < Ndimensions; i++)
    {
        for (int j = 0; j < nestedSampleOfParameters.cols(); j++)
        {
            nestedSampleOfParameters(i,j) = uniform(engine)*(maxima(i)-minima(i)) + minima(i);
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
//      nestedSampleOfParameters: one-dimensional Eigen Array containing the set of 
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

void UniformPrior::drawWithConstraint(RefArrayXd nestedSampleOfParameters, Likelihood &likelihood)
{
    double logLikelihood;
    double logLikelihoodConstraint = likelihood.logValue(nestedSampleOfParameters);


    // Uniform sampling to find new parameter with logLikelihood > logLikelihoodConstraint
    
    do
    {
        for (int i = 0; i < Ndimensions; i++)
            {
                nestedSampleOfParameters(i) = uniform(engine)*(maxima(i) - minima(i)) + minima(i);
            }
    
        logLikelihood = likelihood.logValue(nestedSampleOfParameters);
    }
    while (logLikelihood <= logLikelihoodConstraint);
    

}













// UniformPrior::pointIsRejected()
//
// PURPUSE:
//      Evaluates whether input point coordinates satisfy prior conditions.
//
// INPUT:
//      drawnSampleOfParameters: an Eigen Array of size Ndimensions
//      containing a sample of coordinates for one object to be verified.
//
// OUTPUT:
//      A bool variable declaring whether the point has to be rejected (true)
//      or accepted (false)
//

bool UniformPrior::pointIsRejected(RefArrayXXd drawnSampleOfParameters)
{
    assert (drawnSampleOfParameters.cols() == 1);
    assert (drawnSampleOfParameters.rows() == Ndimensions);


    bool pointIsRejected = false;           // Start without rejection

    for (int i = 0; i < Ndimensions; i++)
    {
        // If at least in one dimension the coordinate is not verifying boundaries, 
        // reject the point and end function

        if ((drawnSampleOfParameters(i,0) > maxima(i)) || (drawnSampleOfParameters(i,0) < minima(i)))
        {
            pointIsRejected = true;
            return pointIsRejected;
        }
    }

    return pointIsRejected;
}
