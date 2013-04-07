// Derived class for identifying intersection between
// different hyper elliposids.
// Created by Enrico Corsaro @ IvS - 27 March 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "HyperEllipsoidIntersector.h"
// Implementation contained in "HyperEllipsoidIntersector.cpp"

#ifndef HYPERELLIPSOIDINTERSECTOR_H
#define HYPERELLIPSOIDINTERSECTOR_H

#include "HyperQuadricIntersector.h"

using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class HyperEllipsoidIntersector : public HyperQuadricIntersector
{
    public:

        HyperEllipsoidIntersector();
        ~HyperEllipsoidIntersector();

        virtual bool intersection(const RefArrayXXd covarianceMatrix1, const RefArrayXd centerCoordinates1, 
                                  const RefArrayXXd covarianceMatrix2, const RefArrayXd centerCoordinates2);
        ArrayXi findNonOverlappingEllipsoids(const int Nclusters, const RefArrayXd allEnlargedEigenvalues, 
                                             const RefArrayXXd allEigenvectorsMatrix, const RefArrayXd allCentersCoordinates);
        bool checkPointForOverlap(const RefArrayXd enlargedEigenValues, const RefArrayXXd eigenVectorsMatrix,
                                  const RefArrayXd centerCoordinates, const RefArrayXd pointCoordinates);

    protected:


    private:

};

#endif
