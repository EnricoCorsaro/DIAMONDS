//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYMETRIC_H
#define DIAMONDS_PYMETRIC_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

using namespace std;

template <class MetricBase = Metric> class PyMetric : public MetricBase
{
public:
    using MetricBase::MetricBase;
    double distance(RefArrayXd point1, RefArrayXd point2) override
    {
        PYBIND11_OVERLOAD_PURE(
                double,
                MetricBase,
                distance,
                point1,
                point2
                );
    }
};

#endif //DIAMONDS_PYMETRIC_H
