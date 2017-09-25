//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYNORMALPRIOR_H
#define DIAMONDS_PYNORMALPRIOR_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "NormalPrior.h"

class PyNormalPrior : public NormalPrior
{
public:
    using NormalPrior::NormalPrior;

    double logDensity(RefArrayXd const x, const bool includeConstantTerm = false) override
    {
        PYBIND11_OVERLOAD(
                double,
                NormalPrior,
                logDensity,
                x,
                includeConstantTerm
        );
    }

    bool drawnPointIsAccepted(RefArrayXd const drawnPoint) override
    {
        PYBIND11_OVERLOAD(
                bool,
                NormalPrior,
                drawnPointIsAccepted,
                drawnPoint
        );
    }

    void draw(RefArrayXXd drawnSample) override
    {
        PYBIND11_OVERLOAD(
                void,
                NormalPrior,
                draw,
                drawnSample
        );
    }

    void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood) override
    {
        PYBIND11_OVERLOAD(
                void,
                NormalPrior,
                drawWithConstraint,
                drawnPoint,
                likelihood
        );
    }

    void writeHyperParametersToFile(string fullPath) override
    {
        PYBIND11_OVERLOAD(
                void,
                NormalPrior,
                writeHyperParametersToFile,
                fullPath
        );
    }
};


#endif //DIAMONDS_PYNORMALPRIOR_H
