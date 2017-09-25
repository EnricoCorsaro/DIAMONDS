//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYSUPERGAUSSIANPRIOR_H
#define DIAMONDS_PYSUPERGAUSSIANPRIOR_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "SuperGaussianPrior.h"

class PySuperGaussianPrior : public SuperGaussianPrior
{
public:
    using SuperGaussianPrior::SuperGaussianPrior;

    double logDensity(RefArrayXd const x, const bool includeConstantTerm = false) override
    {
        PYBIND11_OVERLOAD(
                double,
                SuperGaussianPrior,
                logDensity,
                x,
                includeConstantTerm
        );
    }

    bool drawnPointIsAccepted(RefArrayXd const drawnPoint) override
    {
        PYBIND11_OVERLOAD(
                bool,
                SuperGaussianPrior,
                drawnPointIsAccepted,
                drawnPoint
        );
    }

    void draw(RefArrayXXd drawnSample) override
    {
        PYBIND11_OVERLOAD(
                void,
                SuperGaussianPrior,
                draw,
                drawnSample
        );
    }

    void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood) override
    {
        PYBIND11_OVERLOAD(
                void,
                SuperGaussianPrior,
                drawWithConstraint,
                drawnPoint,
                likelihood
        );
    }

    void writeHyperParametersToFile(string fullPath) override
    {
        PYBIND11_OVERLOAD(
                void,
                SuperGaussianPrior,
                writeHyperParametersToFile,
                fullPath
        );
    }
};


#endif //DIAMONDS_PYSUPERGAUSSIANPRIOR_H
