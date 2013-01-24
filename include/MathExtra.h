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
    
        vector<double> y;
    
        void lorentzProfile(vector<double> &x, double x0, double gamma, double amp = 1)
        {
            unsigned long xsize = x.size();
            y.resize(xsize);

            for (int i=0; i < x.size(); i++)
            {
                y.at(i) = amp/(1 + pow( ((x.at(i)-x0)/gamma), 2) );
            }

            return;
        }
    
    
        void gaussProfile(vector<double> &x, double x0, double sigma, double amp = 1)
        {
            double fac;
            unsigned long xsize = x.size();
            y.resize(xsize);

            fac = 1./(sqrt(2.*PI) * sigma);

            for (int i=0; i < x.size(); i++)
            {
                y.at(i) = fac * exp (-1.*pow( (x.at(i) - x0), 2) / (2.* pow(sigma, 2)));
            }
            return;
        }


        double gaussLikelihood(vector<double> &x_obs, vector<double> &x_theor, vector<double> &sigma)
        {
            if ((x_obs.size() != x_theor.size()) && (x_obs.size() != sigma.size()))
            {
                cout << "Array dimensions do not match. Quitting program." << endl;
                exit(1);
            }
        
            double fac;
            double likel;
            unsigned long xsize = x_obs.size();
            y.resize(xsize);
            vector<double> delta(xsize);

            for (int i=0; i < x_obs.size(); i++)
            {
                fac = 1./(sqrt(2.*PI) * sigma.at(i));
                delta.at(i) = (1./2.) * pow( (x_obs.at(i) - x_theor.at(i)), 2) / 
                (pow(sigma.at(i), 2));
                y.at(i) = fac * exp(-1.*delta.at(i));
            }
    
            likel = product (y);   
            return likel;
        }


        double product(vector<double> &vec)
        {
            unsigned long size = vec.size()     ;
            double prod = 1.    ;
        
            for (int i = 0; i < size; i++)
            {
                prod *= vec.at(i);
            }
            return prod ;
        }
    
    
        double total(vector<double> &vec)
        {
            unsigned long size = vec.size();
            double sum = 0.;
        
            for (int i = 0; i < size; i++)
            {
                sum += vec.at(i);
            }
            return sum  ;
        }
        
    protected:
    
    private:
};
