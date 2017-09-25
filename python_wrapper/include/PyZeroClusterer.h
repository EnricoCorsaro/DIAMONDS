//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYZEROCLUSTERER_H
#define DIAMONDS_PYZEROCLUSTERER_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "ZeroClusterer.h"

class PyZeroClusterer : public ZeroClusterer
{
public:
    using ZeroClusterer::ZeroClusterer;

    int cluster(RefArrayXXd sample, vector<int> &optimalClusterIndices, vector<int> &optimalClusterSizes) override
    {
        PYBIND11_OVERLOAD(
                int,
                ZeroClusterer,
                cluster,
                sample,
                optimalClusterIndices,
                optimalClusterSizes
        );
    }
};

#endif //DIAMONDS_PYZEROCLUSTERER_H
