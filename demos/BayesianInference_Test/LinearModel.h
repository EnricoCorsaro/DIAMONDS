// Derived class for building a linear model
// Created by Enrico Corsaro @ CEA - October 2015
// e-mail: emncorsaro@gmail.com
// Header file "LinearModel.h"
// Implementations contained in "LinearModel.cpp"


#ifndef LINEARMODEL_H
#define LINEARMODEL_H

#include <cstdlib>
#include <Eigen/Core>
#include "Functions.h"
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class LinearModel : public Model
{
    public:
    
        LinearModel(const RefArrayXd covariates);
        ~LinearModel();
        ArrayXd getCovariates();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);

    protected:


    private:
    
}; // END class LinearModel


#endif
