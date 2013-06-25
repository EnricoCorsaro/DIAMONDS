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

        Ellipsoid(RefArrayXXd sampleOfParameters, const double enlargementFactor=1.0);
        ~Ellipsoid();

        void resetEnlargementFactor(const double newEnlargementFactor);
        bool overlapsWith(Ellipsoid ellipsoid);
        bool containsPoint(const RefArrayXd pointCoordinates);
        void drawPoint(RefArrayXd drawnPoint);
        ArrayXd getCenterCoordinates();
        ArrayXd getEigenvalues();
        ArrayXXd getSample();
        ArrayXXd getCovarianceMatrix();
        ArrayXXd getEigenvectors();
        int getSampleSize();
        double getHyperVolume();
        double getEnlargementFactor();


    protected:

        ArrayXd centerCoordinates;
        ArrayXd originalEigenvalues;        // non-enlarged eigenvalues
        ArrayXd enlargedEigenvalues;        // enlarged eigenvalues
        ArrayXXd sample;
        ArrayXXd covarianceMatrix;  
        ArrayXXd eigenvectors;
        int sampleSize;
        double hyperVolume;
        double enlargementFactor;


    private:

        int Ndimensions;
        mt19937 engine;
        uniform_real_distribution<> uniform;
        normal_distribution<> normal;  


};

#endif
