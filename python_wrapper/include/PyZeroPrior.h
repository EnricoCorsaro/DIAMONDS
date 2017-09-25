//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYZEROPRIOR_H
#define DIAMONDS_PYZEROPRIOR_H


#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "ZeroPrior.h"

class PyZeroPrior:public ZeroPrior
{
public:
    using ZeroPrior::ZeroPrior;
    double logDensity(RefArrayXd const x, const bool includeConstantTerm = false) override
    {
        PYBIND11_OVERLOAD(
                double,
                ZeroPrior,
                logDensity,
                x,
                includeConstantTerm
        );
    }

    bool drawnPointIsAccepted(RefArrayXd const drawnPoint) override
    {
        PYBIND11_OVERLOAD(
                bool,
                ZeroPrior,
                drawnPointIsAccepted,
                drawnPoint
        );
    }

    void draw(RefArrayXXd drawnSample) override
    {
        PYBIND11_OVERLOAD(
                void,
                ZeroPrior,
                draw,
                drawnSample
        );
    }

    void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood) override
    {
        PYBIND11_OVERLOAD(
                void,
                ZeroPrior,
                drawWithConstraint,
                drawnPoint,
                likelihood
        );
    }

    void writeHyperParametersToFile(string fullPath) override
    {
        PYBIND11_OVERLOAD(
                void,
                ZeroPrior,
                writeHyperParametersToFile,
                fullPath
        );
    }
};

#endif //DIAMONDS_PYZEROPRIOR_H
