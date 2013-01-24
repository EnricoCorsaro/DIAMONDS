// Class for implementing some extra functions useful for Peak Bagging fitting process
// Created by Enrico Corsaro @ IvS - 22 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "MathExtra.h"
#include "lib.h"
#define PI 3.14159265358979
using namespace std;

class MathExtra
{
    public:
    
        // Public data members accessible from outside
        vector<double> y;
	
        // Public member functions
        // Computes a simple Lorentzian profile given the centroid and the width
        void lorentzProfile ( vector<double> &x, double x0, double gamma, double amp = 1 )
        {
            unsigned long xsize = x.size();

            y.resize(xsize);
	    
            for ( int i=0; i < x.size(); i++ )
            {
                y.at(i) = amp/(1 + pow( ((x.at(i)-x0)/gamma), 2) );
            }

            return ;
        }

        // Computes a simple Gaussian profile given the centroid, the standard deviation and the amplitude
        void gaussProfile ( vector<double> &x, double x0, double sigma, double amp = 1 )
        {
            double fac;
            unsigned long xsize = x.size();

            y.resize(xsize);
            fac = 1./(sqrt(2.*PI) * sigma);

            for ( int i = 0; i < x.size(); i++ )
            {
                y.at(i) = fac * exp (-1.*pow( (x.at(i) - x0), 2) / (2.* pow(sigma, 2)));
            }

            return ;
        }

        // Computes the gaussian likelihood from a set of observations, uncertainties,
        // and theoretical predictions
        double gaussLikelihood ( vector<double> &x_obs, vector<double> &x_theor, vector<double> &sigma )
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

            for ( int i = 0; i < xsize; i++ )
            {
                fac = 1./(sqrt(2.*PI) * sigma.at(i));
                delta.at(i)	= -0.5*pow( (x_obs.at(i) - x_theor.at(i)), 2) /
				(pow(sigma.at(i), 2));
                likelihood.at(i) = fac * exp(delta.at(i));
             }
	
            likel = product (likelihood);

            return likel;
        }
	
        // Computes the logarithmic gaussian likelihood from a set of observations, uncertainties,
        // and theoretical predictions
        double logGaussLikelihood ( vector<double> &x_obs, vector<double> &x_theor, vector<double> &sigma )
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

            likel = total (lambda);

            return likel;
        }

        // Computes the product of the elements contained in a vector of doubles
        double product ( vector<double> &vec )
        {
            unsigned long size = vec.size();
            double prod = 1.;
		
            for (unsigned long i = 0; i < size; i++)
            {
                prod *= vec.at(i);
            }

            return prod;
        }

        // Computes the sum of the elements contained in a vector of doubles
        double total ( vector<double> &vec )
        {
            unsigned long size = vec.size();
            double sum = 0.;
		
            for (unsigned long i = 0; i < size; i++)
            {
                sum += vec.at(i);
            }

            return sum;
        }

        // Computes a logaritmic summation of exponentials in order to avoid overflow errors
        double logExpSum ( double &x, double &y )
        {
            if (x >= y )
            {
                return x + log(1. + exp(y - x));
            }
            else
            {
                return y + log(1. + exp(x - y));
            }
        }
} ;
