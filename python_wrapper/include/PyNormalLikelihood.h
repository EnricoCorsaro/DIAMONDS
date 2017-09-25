//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYNORMALLIKELIHOOD_H
#define DIAMONDS_PYNORMALLIKELIHOOD_H
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "NormalLikelihood.h"

class PyNormalLikelihood : public NormalLikelihood
{
public:
    using NormalLikelihood::NormalLikelihood;

    double logValue(RefArrayXd const modelParameters) override
    {
        PYBIND11_OVERLOAD(
                double,
                NormalLikelihood,
                logValue,
                modelParameters
        );
    }
};

#endif //DIAMONDS_PYNORMALLIKELIHOOD_H
