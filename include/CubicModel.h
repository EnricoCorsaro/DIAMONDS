// Derived class for building a linear model
// Created by Enrico Corsaro @ OACT - August 2017
// e-mail: emncorsaro@gmail.com
// Header file "CubicModel.h"
// Implementations contained in "CubicModel.cpp"


#ifndef LINEARMODEL_H
#define LINEARMODEL_H

#include <cstdlib>
#include <Eigen/Core>
#include "Functions.h"
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class CubicModel : public Model
{
    public:
    
        CubicModel(const RefArrayXd covariates);
        ~CubicModel();
        ArrayXd getCovariates();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);

    protected:


    private:
    
}; // END class CubicModel


#endif
