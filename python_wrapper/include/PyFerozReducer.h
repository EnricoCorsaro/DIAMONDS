//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYFEROZREDUCER_H
#define DIAMONDS_PYFEROZREDUCER_H
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "FerozReducer.h"

class PyFerozReducer : public FerozReducer
{
public:
    using FerozReducer::FerozReducer;
    int updateNlivePoints() override
    {
        PYBIND11_OVERLOAD(
                int,
                FerozReducer,
                updateNlivePoints
        );
    }
};

#endif //DIAMONDS_PYFEROZREDUCER_H
