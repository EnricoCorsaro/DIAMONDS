//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYPRIOR_H
#define DIAMONDS_PYPRIOR_H
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "Prior.h"

class PyPrior:public Prior
{
public:
    using Prior::Prior;
    double logDensity(RefArrayXd const x, const bool includeConstantTerm = false) override
    {
        PYBIND11_OVERLOAD_PURE(
                double,
                Prior,
                logDensity,
                x,
                includeConstantTerm
        );
    }

    bool drawnPointIsAccepted(RefArrayXd const drawnPoint) override
    {
        PYBIND11_OVERLOAD_PURE(
                bool,
                Prior,
                drawnPointIsAccepted,
                drawnPoint
        );
    }

    void draw(RefArrayXXd drawnSample) override
    {
        PYBIND11_OVERLOAD_PURE(
                void,
                Prior,
                draw,
                drawnSample
        );
    }

    void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood) override
    {
        PYBIND11_OVERLOAD_PURE(
                void,
                Prior,
                drawWithConstraint,
                drawnPoint,
                likelihood
        );
    }

    void writeHyperParametersToFile(string fullPath) override
    {
        PYBIND11_OVERLOAD_PURE(
                void,
                Prior,
                writeHyperParametersToFile,
                fullPath
        )
    }
};

#endif //DIAMONDS_PYPRIOR_H
