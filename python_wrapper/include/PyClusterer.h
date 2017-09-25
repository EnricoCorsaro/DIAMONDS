//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYCLUSTERER_H
#define DIAMONDS_PYCLUSTERER_H
#include "Clusterer.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;

class PyClusterer : public Clusterer
{
public:
    using Clusterer::Clusterer;
    int cluster(RefArrayXXd sample, vector<int> &optimalClusterIndices, vector<int> &optimalClusterSizes) override {
        PYBIND11_OVERLOAD_PURE(
        int,
        Clusterer,
        cluster,
        sample,
        optimalClusterIndices,
        optimalClusterSizes
        );
    }
};

#endif //DIAMONDS_PYCLUSTERER_H
