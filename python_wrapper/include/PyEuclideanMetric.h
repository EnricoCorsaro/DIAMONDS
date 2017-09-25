//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYEUCLIDEANMETRIC_H
#define DIAMONDS_PYEUCLIDEANMETRIC_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "EuclideanMetric.h"

using namespace std;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;

class PyEuclideanMetric : public EuclideanMetric
{
public:
    using EuclideanMetric::EuclideanMetric;
    double distance(RefArrayXd point1,RefArrayXd point2) override
    {
        PYBIND11_OVERLOAD(
                double,
                EuclideanMetric,
                distance,
                point1,
                point2
        );
    }

};

#endif //DIAMONDS_PYEUCLIDEANMETRIC_H
