// Configuration code for setting up nested sampling analysis
// Created by Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Source code file "configure.cpp"
#include "MathExtra.h"
#include "FileProcess.h"

int main()
{
	double xmax,xmin	;	// Maximum and minimum parameter values
	vector<double> uniform	;	// vector of uniform numbers extracted in the range (0,1)
	vector<double> x	;	// vector of x numbers extracted from uniform
	double prior	;	// 1-dimensional flat prior constant value
	int n		; 	// Number of objects to be extracted
	MathExtra	math	;	// Object providing extra math functions

	n	= 100;
	xmax	= 20;
	xmin	= 0;
	uniform.resize(n);
	x.resize(n);
	prior	= pow((xmax - xmin),-1.);
	srand(time(0));
	
	for (int i=0; i < n; i++)
	{
		uniform.at(i) 	= rand()/(RAND_MAX + 1.)	;
		x.at(i) 	= uniform.at(i)/prior + xmin	;
		//cout << x.at(i) << endl;
	}

	// Example of a Gaussian likelihood
	double x0 = 10;
	double sigma = 2.5;

	math.gaussProfile ( x, x0, sigma );	// output stored in math.y

	
}
