//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYNESTEDSAMPLER_H
#define DIAMONDS_PYNESTEDSAMPLER_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "NestedSampler.h"

class PyNestedSampler : public NestedSampler
{
public:
    using NestedSampler::NestedSampler;

    bool drawWithConstraint(const RefArrayXXd totalSample, const unsigned int Nclusters, const vector<int> &clusterIndices,
                                    const vector<int> &clusterSizes, RefArrayXd drawnPoint,
                                    double &logLikelihoodOfDrawnPoint, const int maxNdrawAttempts) override
    {
        PYBIND11_OVERLOAD_PURE(
        bool,
        NestedSampler,
        drawWithConstraint,
        totalSample,
        Nclusters,
        clusterIndices,
        clusterSizes,
        drawnPoint,
        logLikelihoodOfDrawnPoint,
        maxNdrawAttempts
        );

    }

    bool verifySamplerStatus() override
    {
        PYBIND11_OVERLOAD_PURE(
                bool,
                NestedSampler,
                verifySamplerStatus
        );
    }

};

class NestedSamplerPublicist : public NestedSampler {
public:
    using NestedSampler::verifySamplerStatus; // inherited with different access modifier
};


#endif //DIAMONDS_PYNESTEDSAMPLER_H
