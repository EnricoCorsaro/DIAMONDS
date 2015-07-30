// Derived class for zero model objects. 
// This model returns zero for any set of model parameters.
// Created by Enrico Corsaro @ IvS - 9 May 2014
// e-mail: emncorsaro@gmail.com
// Header file "ZeroModel.h"
// Implementations contained in "ZeroModel.cpp"


#ifndef ZEROMODEL_H
#define ZEROMODEL_H

#include <iostream>
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class ZeroModel : public Model
{
    public:
    
        ZeroModel(const RefArrayXd covariates);
        ~ZeroModel();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters) override;
       
    protected:


    private:
    
};


#endif
