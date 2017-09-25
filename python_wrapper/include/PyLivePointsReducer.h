//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYLIVEPOINTSREDUCER_H
#define DIAMONDS_PYLIVEPOINTSREDUCER_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "LivePointsReducer.h"

class PyLivePointsReducer : public LivePointsReducer
{
public:
    using LivePointsReducer::LivePointsReducer;
    int updateNlivePoints() override
    {
        PYBIND11_OVERLOAD_PURE(
                int,
                LivePointsReducer,
                updateNlivePoints
        );
    }
};

#endif //DIAMONDS_PYLIVEPOINTSREDUCER_H
