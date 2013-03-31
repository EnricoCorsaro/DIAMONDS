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

        virtual bool intersection(RefArrayXXd covarianceMatrix1, RefArrayXXd covarianceMatrix2, 
                                    RefArrayXd centerCoordinates1, RefArrayXd centerCoordinates2);


    protected:


    private:

};

#endif
