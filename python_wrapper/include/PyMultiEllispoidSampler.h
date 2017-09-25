//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYMULTIELLISPOIDSAMPLER_H
#define DIAMONDS_PYMULTIELLISPOIDSAMPLER_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "MultiEllipsoidSampler.h"

class PyMultiEllipsoidSampler : public MultiEllipsoidSampler
{
public:
    using MultiEllipsoidSampler::MultiEllipsoidSampler;
    bool drawWithConstraint(const RefArrayXXd totalSample, const unsigned int Nclusters, const vector<int> &clusterIndices,
                                    const vector<int> &clusterSizes, RefArrayXd drawnPoint,
                                    double &logLikelihoodOfDrawnPoint, const int maxNdrawAttempts) override
    {
        PYBIND11_OVERLOAD(
                bool,
                MultiEllipsoidSampler,
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
                MultiEllipsoidSampler,
                verifySamplerStatus
        );
    }
};

#endif //DIAMONDS_PYMULTIELLISPOIDSAMPLER_H
