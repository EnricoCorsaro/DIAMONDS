//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYMEANNORMALLIKELIHOOD_H
#define DIAMONDS_PYMEANNORMALLIKELIHOOD_H
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "MeanNormalLikelihood.h"

class PyMeanNormalLikelihood : public MeanNormalLikelihood
{
public:
    using MeanNormalLikelihood::MeanNormalLikelihood;
    double logValue(RefArrayXd const modelParameters) override
    {
        PYBIND11_OVERLOAD(
                double,
                MeanNormalLikelihood,
                logValue,
                modelParameters
        );
    }
};


#endif //DIAMONDS_PYMEANNORMALLIKELIHOOD_H
