//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYMETRIC_H
#define DIAMONDS_PYMETRIC_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

using namespace std;

class PyMetric : public Metric
{
public:
    using Metric::Metric;
    double distance(RefArrayXd point1, RefArrayXd point2) override
    {
        PYBIND11_OVERLOAD_PURE(
                double,
                Metric,
                distance,
                point1,
                point2
                );
    }
};

#endif //DIAMONDS_PYMETRIC_H
