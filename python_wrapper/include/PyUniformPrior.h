//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYUNIFORMPRIOR_H
#define DIAMONDS_PYUNIFORMPRIOR_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "UniformPrior.h"

class PyUniformPrior : public UniformPrior
{
public:
    using UniformPrior::UniformPrior;

    double logDensity(RefArrayXd const x, const bool includeConstantTerm = false) override
    {
        PYBIND11_OVERLOAD(
                double,
                UniformPrior,
                logDensity,
                x,
                includeConstantTerm
        );
    }

    bool drawnPointIsAccepted(RefArrayXd const drawnPoint) override
    {
        PYBIND11_OVERLOAD(
                bool,
                UniformPrior,
                drawnPointIsAccepted,
                drawnPoint
        );
    }

    void draw(RefArrayXXd drawnSample) override
    {
        PYBIND11_OVERLOAD(
                void,
                UniformPrior,
                draw,
                drawnSample
        );
    }

    void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood) override
    {
        PYBIND11_OVERLOAD(
                void,
                UniformPrior,
                drawWithConstraint,
                drawnPoint,
                likelihood
        );
    }

    void writeHyperParametersToFile(string fullPath) override
    {
        PYBIND11_OVERLOAD(
                void,
                UniformPrior,
                writeHyperParametersToFile,
                fullPath
        );
    }
};

#endif //DIAMONDS_PYUNIFORMPRIOR_H
