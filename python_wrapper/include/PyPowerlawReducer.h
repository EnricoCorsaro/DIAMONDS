//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYPOWERLAWREDUCER_H
#define DIAMONDS_PYPOWERLAWREDUCER_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "PowerlawReducer.h"

class PyPowerlawReducer : public PowerlawReducer
{
public:
    using PowerlawReducer::PowerlawReducer;

    int updateNlivePoints() override
    {
        PYBIND11_OVERLOAD(
                int,
                PowerlawReducer,
                updateNlivePoints
        );
    }

};

#endif //DIAMONDS_PYPOWERLAWREDUCER_H
