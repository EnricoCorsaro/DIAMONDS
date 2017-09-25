//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYZEROSAMPLER_H
#define DIAMONDS_PYZEROSAMPLER_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "ZeroSampler.h"

class PyZeroSampler: public ZeroSampler
{
public:
    using ZeroSampler::ZeroSampler;
    bool drawWithConstraint(const RefArrayXXd totalSample, const unsigned int Nclusters, const vector<int> &clusterIndices,
                                    const vector<int> &clusterSizes, RefArrayXd drawnPoint,
                                    double &logLikelihoodOfDrawnPoint, const int maxNdrawAttempts) override
    {
        PYBIND11_OVERLOAD(
                bool,
                ZeroSampler,
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
        PYBIND11_OVERLOAD(
                bool,
                ZeroSampler,
                verifySamplerStatus
        );
    }
};

#endif //DIAMONDS_PYZEROSAMPLER_H
