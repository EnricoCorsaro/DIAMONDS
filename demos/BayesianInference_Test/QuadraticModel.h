// Derived class for building a linear model
// Created by Enrico Corsaro @ OACT - August 2017
// e-mail: emncorsaro@gmail.com
// Header file "QuadraticModel.h"
// Implementations contained in "QuadraticModel.cpp"


#ifndef LINEARMODEL_H
#define LINEARMODEL_H

#include <cstdlib>
#include <Eigen/Core>
#include "Functions.h"
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class QuadraticModel : public Model
{
    public:
    
        QuadraticModel(const RefArrayXd covariates);
        ~QuadraticModel();
        ArrayXd getCovariates();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);

    protected:


    private:
    
}; // END class QuadraticModel


#endif
