// Class for creating an ellipsoid object to be
// used within the sampler class.
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 11 April 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "Ellipsoid.h"
// Implementation contained in "Ellipsoid.cpp"

#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include <iostream>
#include <iomanip>
#include <cfloat>
#include <random>
#include <ctime>
#include <cmath>
#include <vector>
#include <cassert>
#include <Eigen/Dense>
#include "Functions.h"

using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class Ellipsoid
{

    public:

        Ellipsoid(RefArrayXXd sampleOfParameters, const int index);
        ~Ellipsoid();

        void build(const double enlargementFactor);
        ArrayXd getCenterCoordinates();
        ArrayXd getEigenvalues();
        ArrayXXd getSampleOfParameters();
        ArrayXXd getCovarianceMatrix();
        ArrayXXd getEigenvectorsMatrix();
        int getIndex();
        int getNobjects();
        double getHyperVolume();
        double getEnlargementFactor();


    protected:

        ArrayXd centerCoordinates;
        ArrayXd eigenvalues;
        ArrayXXd sampleOfParameters;
        ArrayXXd covarianceMatrix;
        ArrayXXd eigenvectorsMatrix;
        int Nobjects;
        int index;
        double hyperVolume;
        double enlargementFactor;


    private:

        int Ndimensions;

};

#endif
