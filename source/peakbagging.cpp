// Main code for peak bagging by means of nested sampling analysis
// Created by Enrico Corsaro @ IvS - 24 January 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Source code file "peakbagging.cpp"
#include "MathExtra.h"
#include "FileProcess.h"

int main()
{
	double logwidth;	// ln(width in prior mass)
	double H = 0.;		// Information achieved from prior to posterior, initially  0
	double logZ = 0.;	// ln(Evidence Z), initially 0
	double logZnew;		// Updated logZ
	int nest;           // Nested sampling iteration count

	// Outermost interval of prior mass: width = 1 means total prior mass
	logwidth = log(1. - exp(-1. / n));
}
