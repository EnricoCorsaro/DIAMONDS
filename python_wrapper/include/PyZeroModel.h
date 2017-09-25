//
// Created by mamu on 9/22/17.
//

#ifndef DIAMONDS_PYZEROMODEL_H
#define DIAMONDS_PYZEROMODEL_H

#include "ZeroModel.h"

class PyZeroModel : public ZeroModel
{
public:
    using ZeroModel::ZeroModel;
    void predict(RefArrayXd predictions, const RefArrayXd modelParameters) override
    {
        PYBIND11_OVERLOAD(
                void,
                ZeroModel,
                predict,
                predictions,
                modelParameters
        );
    }
};

#endif //DIAMONDS_PYZEROMODEL_H
