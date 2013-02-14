#include "MathExtra.h"



// MathExtra::lorentzProfile()
//
// PURPOSE: 
//      Computes a simple Lorentzian profile given the centroid and the width.
//      Saves the dependent variable values into a vector y accessible as
//      public data member.
//
// INPUT:
//      y : vector containing the result
//      x : vector containing independent variable values
//      x0 : centroid of the Lorentzian profile
//      gamma : width of the Lorentzian profile (mode linewidth)
//      amp : maximum height of the Lorentzian profile
//
// OUTPUT:
//      void

void MathExtra::lorentzProfile(RefArrayXd y, const RefArrayXd x, double x0, double gamma, double amp)
{
    
        y = (amp*amp)/((x-x0)*(x-x0) + (gamma/2.)*(gamma/2.));

    return;
} // END MathExtra::lorentzProfile()






// MathExtra::logGaussProfile()
//
// PURPOSE: 
//      Computes the logarithm of a Gaussian profile given the centroid, 
//      the standard deviation and the amplitude, for one given x-value
//
// INPUT:
//      x : independent variable
//      mu : mean value of the Gaussian profile
//      sigma : standard deviation of the Gaussian profile
//      amp : maximum amplitude of the Gaussian profile (default = 1)
//
// OUTPUT:
//      logGaussian of the independent variable

double MathExtra::logGaussProfile(double x, double mu, double sigma, double amp)
{
    const double prefactor = log(amp) - 0.5 * log(2*MathExtra::PI) - log(sigma);
    return prefactor - 0.5 * (x - mu) * (x - mu) / (sigma * sigma);
} 




// MathExtra::logGaussProfile()
//
// PURPOSE: 
//      Computes the logarithm of a Gaussian profile given the centroid, 
//      the standard deviation and the amplitude for a set of values.
//
// INPUT:
//      y : Eigen array containing the result
//      x : Eigen array containing independent variable values
//      mu : mean value of the Gaussian profile
//      sigma : standard deviation of the Gaussian profile
//      amp : maximum amplitude of the Gaussian profile (default = 1)
//
// OUTPUT:
//      void
//

void MathExtra::logGaussProfile(RefArrayXd y, const RefArrayXd x, const double mu, const double sigma, const double amp)
{
    const double prefactor = log(amp) - 0.5 * log(2*MathExtra::PI) - log(sigma);
    y = prefactor - 0.5 * (x - mu) * (x - mu) / (sigma * sigma);
    
    return;
} 





// MathExtra::gaussLikelihood()
//
// PURPOSE: 
//      Computes the gaussian likelihood from a set of observations,
//      uncertainties, and theoretical predictions.
//
// INPUT:
//      x_obs = vector containing the oberved values
//      x_theor = vector containing the corresponding predicted values of x
//      sigma = vector containing the uncertanties on the observed values
//
// OUTPUT:
//      The resulting value of the Gaussian likelihood
//

double MathExtra::gaussLikelihood(const RefArrayXd x_obs, const RefArrayXd x_theor, const RefArrayXd sigma)
{
    if ((x_obs.size() != x_theor.size()) && (x_obs.size() != sigma.size()))
    {
        cout << "Array dimensions do not match. Quitting program." << endl;
        exit(1);
    }
    
    double likel;
    ArrayXd fac;
    ArrayXd delta;
    ArrayXd likelihood;
    
    fac = 1./(sqrt(2.*MathExtra::PI) * sigma);
    delta = -0.5*((x_obs - x_theor)*(x_obs - x_theor)) / (sigma*sigma);
    likelihood = fac * delta.exp();
    likel = product(likelihood);
    
    return likel;
} // END MathExtra::gaussianLikelihood()




// MathExtra::logGaussLikelihood() 
//
// PURPOSE: 
//      Computes the log-Gaussian likelihood from a set of observations,
//      uncertainties, and theoretical predictions.
//
// INPUT:
//      x_obs : vector containing the oberved values
//      x_theor : vector containing the corresponding predicted values of x
//      sigma : vector containing the uncertanties on the observed values
//
// OUTPUT:
//      The resulting value of the log-Gaussian likelihood
//

double MathExtra::logGaussLikelihood(const RefArrayXd x_obs, const RefArrayXd x_theor, const RefArrayXd sigma)
{
    if ((x_obs.size() != x_theor.size()) && (x_obs.size() != sigma.size()))
    {
        cout << "Array dimensions do not match. Quitting program." << endl;
        exit(1);
    }
    
    ArrayXd delta;
    ArrayXd lambda0;
    ArrayXd lambda;
    
    delta = ((x_obs - x_theor)*(x_obs - x_theor)) / (sigma*sigma);
    lambda0 = -1.*log(sqrt(2.*PI) * sigma);
    lambda = lambda0 -0.5*delta;
    
    return sum(lambda);

} // END MathExtra::logGaussLikelihood()




// MathExtra::product()
// PURPOSE: 
//      Computes the product of the elements contained in a vector of doubles of class vector.
// INPUT:
//      vec : vector of values to be multiplied
// OUTPUT:
//      The productoria of the vector elements
//

inline double MathExtra::product(const vector<double> &vec)
{
    return accumulate(vec.begin(), vec.end(), 1.0, multiplies<double>());
} // END MathExtra::product()





// MathExtra::product()
// PURPOSE: 
//      Computes the product of the elements contained in a vector of doubles of class Eigen Array.
// INPUT:
//      vec : vector of values to be multiplied
// OUTPUT:
//      The productoria of the vector elements
//
// REMARK: 
//      Is there any member function of this class allowing the sum of the elements? I am not finding it.

double MathExtra::product(const RefArrayXd vec)
{
    double prod = 1.;
    for (ptrdiff_t i = 0; i < vec.size(); i++)
    {
        prod *= vec(i);
    }
    return prod;
} // END MathExtra::product()




// MathExtra::sum()
// PURPOSE: 
//      Computes the sum of the elements contained in a vector of doubles of vector class
// INPUT:
//      vec = vector of values to be added
// OUTPUT:
//      The summation of the vector elements
//

inline double MathExtra::sum(const vector<double> &vec)
{
    return accumulate(vec.begin(), vec.end(), 0.0);
} // END MathExtra::sum()




// MathExtra::sum()
// PURPOSE: 
//      Computes the sum of the elements contained in a vector of doubles of Eigen Array class
// INPUT:
//      vec = vector of values to be added
// OUTPUT:
//      The summation of the vector elements
//
// REMARK: 
//      Is there any member function of this class allowing the sum of the elements? I am not finding it.

double MathExtra::sum(const RefArrayXd vec)
{
    double total = 0.;
    for (ptrdiff_t i = 0; i < vec.size(); i++)
    {
        total += vec(i);
    }
    return total;
} // END MathExtra::sum()




// MathExtra::logExpSum()
// PURPOSE: 
//      Computes a logaritmic summation of exponentials in order to avoid overflow errors
// INPUT: 
//      x = a first variable to be added
//      y = a second variable to be added
// OUTPUT: 
//        The logarithmic summation of the exponentials log(exp(x)+exp(y))
//

double MathExtra::logExpSum(double x, double y)
{
    return (x >= y ? x+log(1.+exp(y-x)) : y + log(1.+exp(x-y)));
} // END MathExtra::logExpSum()

