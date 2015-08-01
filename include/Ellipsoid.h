// Class for creating an ellipsoid object to be
// used within the sampler class.
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 11 April 2013
// e-mail: emncorsaro@gmail.com
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

        Ellipsoid(RefArrayXXd sampleOfParameters, const double enlargementFraction = 0.0);              // Default no enlargement
        ~Ellipsoid();

        void resetEnlargementFraction(const double newEnlargementFraction);
        bool overlapsWith(Ellipsoid ellipsoid, bool &ellipsoidMatrixDecompositionIsSuccessful);
        bool containsPoint(const RefArrayXd pointCoordinates);
        void drawPoint(RefArrayXd drawnPoint);
        ArrayXd getCenterCoordinates();
        ArrayXd getEigenvalues();
        ArrayXXd getSample();
        ArrayXXd getCovarianceMatrix();
        ArrayXXd getEigenvectors();
        int getSampleSize();
        double getHyperVolume();
        double getEnlargementFraction();


    protected:

        ArrayXd centerCoordinates;
        ArrayXd originalEigenvalues;        // non-enlarged eigenvalues
        ArrayXd enlargedEigenvalues;        // enlarged eigenvalues
        ArrayXXd sample;
        ArrayXXd covarianceMatrix;  
        ArrayXXd eigenvectors;
        int sampleSize;
        double hyperVolume;
        double enlargementFraction;


    private:

        int Ndimensions;
        mt19937 engine;
        uniform_real_distribution<> uniform;
        normal_distribution<> normal;  


};

#endif
