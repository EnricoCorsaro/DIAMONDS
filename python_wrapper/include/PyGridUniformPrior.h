//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYGRIDUNIFORMPRIOR_H
#define DIAMONDS_PYGRIDUNIFORMPRIOR_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "GridUniformPrior.h"

using namespace std;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;

class PyGridUniformPrior : public GridUniformPrior
{
public:
    using GridUniformPrior::GridUniformPrior;
    double logDensity(RefArrayXd const x, const bool includeConstantTerm = false) override
    {
        PYBIND11_OVERLOAD(
                double,
                GridUniformPrior,
                logDensity,
                x,
                includeConstantTerm
        );
    }
    bool drawnPointIsAccepted(RefArrayXd const drawnPoint) override
    {
        PYBIND11_OVERLOAD(
                bool,
                GridUniformPrior,
                drawnPointIsAccepted,
                drawnPoint
        );
    }
    void draw(RefArrayXXd drawnSample) override {
        PYBIND11_OVERLOAD(
                void,
                GridUniformPrior,
                draw,
                drawnSample
        );
    }
    void drawWithConstraint(RefArrayXd drawnPoint, Likelihood &likelihood) override
    {
        PYBIND11_OVERLOAD(
                void,
                GridUniformPrior,
                drawWithConstraint,
                drawnPoint,
                likelihood
        );
    }
    void writeHyperParametersToFile(string fullPath) override
    {
        PYBIND11_OVERLOAD(
                void,
                GridUniformPrior,
                writeHyperParametersToFile,
                fullPath
        );
    }
};

#endif //DIAMONDS_PYGRIDUNIFORMPRIOR_H
