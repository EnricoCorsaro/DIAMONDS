// Base class for identifying intersection between
// different hyper quadric surfaces.
// Created by Enrico Corsaro @ IvS - 27 March 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "HyperQuadricIntersector.h"
// Implementation contained in "HyperQuadricIntersector.cpp"

#ifndef HYPERQUADRICINTERSECTOR_H
#define HYPERQUADRICINTERSECTOR_H

#include <Eigen/Dense>
#include <cassert>
#include <iostream>

using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class HyperQuadricIntersector
{
    public:

        HyperQuadricIntersector();
        ~HyperQuadricIntersector();

        virtual bool intersection(RefArrayXXd quadricMatrix1, RefArrayXXd quadricMatrix2, 
                                  RefArrayXd centerCoordinates1, RefArrayXd centerCoordinates2) = 0;


    protected:


    private:

};

#endif
