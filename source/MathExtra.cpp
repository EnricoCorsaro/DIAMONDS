
#include "MathExtra.h"

// MathExtra::
//
// PURPOSE: 
//
// INPUT:
//
// OUTPUT:
//




// MathExtra::lorentzProfile()
//
// PURPOSE: Computes a simple Lorentzian profile given the centroid and the width
//
// INPUT:
//
// OUTPUT:
//

void MathExtra::lorentzProfile(vector<double> &x, double x0, double gamma, double amp)
{
    unsigned long xsize = x.size();
    
    y.resize(xsize);
    
    for (unsigned int i=0; i < x.size(); i++)
    {
        y.at(i) = amp/(1 + pow( ((x.at(i)-x0)/gamma), 2) );
    }
    
    return;
}






// MathExtra::gaussProfile()
//
// PURPOSE: Computes a simple Gaussian profile given the centroid, the standard deviation and the amplitude
//
// INPUT:
//
// OUTPUT:
//

void MathExtra::gaussProfile(vector<double> &x, double x0, double sigma, double amp)
{
    double fac;
    unsigned long xsize = x.size();
    
    y.resize(xsize);
    fac = 1./(sqrt(2.*PI) * sigma);
    
    for (unsigned int i = 0; i < x.size(); i++)
    {
        y.at(i) = fac * exp (-1.*pow( (x.at(i) - x0), 2) / (2.* pow(sigma, 2)));
    }
    
    return;
}





// MathExtra::gaussLikelihood()
//
// PURPOSE: Computes the gaussian likelihood from a set of observations,
//          uncertainties, and theoretical predictions
//
// INPUT:
//
// OUTPUT:
//

double MathExtra::gaussLikelihood(vector<double> &x_obs, vector<double> &x_theor, vector<double> &sigma)
{
    if ((x_obs.size() != x_theor.size()) && (x_obs.size() != sigma.size()))
    {
        cout << "Array dimensions do not match. Quitting program." << endl;
        exit(1);
    }
    
    double fac;
    double likel;
    unsigned long xsize = x_obs.size();
    vector<double> delta(xsize);
    vector<double> likelihood(xsize);
    
    for (unsigned int i = 0; i < xsize; i++)
    {
        fac = 1./(sqrt(2.*PI) * sigma.at(i));
        delta.at(i)	= -0.5*pow( (x_obs.at(i) - x_theor.at(i)), 2) /
        (pow(sigma.at(i), 2));
        likelihood.at(i) = fac * exp(delta.at(i));
    }
	
    likel = product (likelihood);
    
    return likel;
}






// MathExtra::logGaussLikelihood() 
//
// PURPOSE: Computes the logarithmic gaussian likelihood from a set of observations,
//          uncertainties, and theoretical predictions
//
// INPUT:
//
// OUTPUT:
//

double MathExtra::logGaussLikelihood(vector<double> &x_obs, vector<double> &x_theor, vector<double> &sigma)
{
    if ((x_obs.size() != x_theor.size()) && (x_obs.size() != sigma.size()))
    {
        cout << "Array dimensions do not match. Quitting program." << endl;
        exit(1);
    }
    
    double likel;
    unsigned long xsize = x_obs.size();
    vector<double> delta(xsize);
    vector<double> lambda0(xsize);
    vector<double> lambda(xsize);
    
    for (unsigned long j = 0; j < xsize; j++)
    {
        delta.at(j)	= pow( ( x_obs.at(j) - x_theor.at(j) ), 2) / ( pow(sigma.at(j), 2) );
        lambda0.at(j) = -1.*log(sqrt(2.*PI) * sigma.at(j));
        lambda.at(j) = lambda0.at(j) -0.5*delta.at(j);
    }
    
    likel = sum(lambda);
    
    return likel;
}





// MathExtra::product()
//
// PURPOSE: Computes the product of the elements contained in a vector of doubles
//
// INPUT:
//
// OUTPUT:
//

inline double MathExtra::product(const vector<double> &vec)
{
    return accumulate(vec.begin(), vec.end(), 1.0, multiplies<double>());
}




// MathExtra::sum()
//
// PURPOSE: Computes the sum of the elements contained in a vector of doubles
//
// INPUT:
//
// OUTPUT:
//

inline double MathExtra::sum(const vector<double> &vec)
{
    return accumulate(vec.begin(), vec.end(), 0.0);
}





// MathExtra::logExpSum()
//
// PURPOSE: Computes a logaritmic summation of exponentials in order to avoid overflow errors
//
// INPUT: x (double)
//        y (double)
//
// OUTPUT: log(exp(x)+exp(y))
//

inline double MathExtra::logExpSum(const double x, const double y)
{
    return (x>=y ? x+log(1.+exp(y-x)) : y + log(1.+exp(x-y)));
}

