//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYKMEANSCLUSTERER_H
#define DIAMONDS_PYKMEANSCLUSTERER_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "KmeansClusterer.h"


using namespace std;

class PyKmeansClusterer : public KmeansClusterer
{
public:
    using KmeansClusterer::KmeansClusterer;
    int cluster(RefArrayXXd sample, vector<int> &optimalClusterIndices, vector<int> &optimalClusterSizes) override
    {
        PYBIND11_OVERLOAD(
                int,
                KmeansClusterer,
                cluster,
                sample,
                optimalClusterIndices,
                optimalClusterSizes
        );
    }
};



#endif //DIAMONDS_PYKMEANSCLUSTERER_H
