//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYLIKELIHOOD_H
#define DIAMONDS_PYLIKELIHOOD_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "Likelihood.h"

class PyLikelihood : public Likelihood
{
public:
    using Likelihood::Likelihood;
    double logValue(RefArrayXd const modelParameters) override
    {
        PYBIND11_OVERLOAD_PURE(
                double,
                Likelihood,
                logValue,
                modelParameters
        );
    }

};

#endif //DIAMONDS_PYLIKELIHOOD_H
