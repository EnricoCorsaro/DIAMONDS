// Derived class for building a super Gaussian model
// Created by Enrico Corsaro @ OACT - February 2017
// e-mail: emncorsaro@gmail.com
// Header file "SuperGaussianModel.h"
// Implementations contained in "SuperGaussianModel.cpp"


#ifndef SUPERGAUSSIANMODEL_H
#define SUPERGAUSSIANMODEL_H

#include <cstdlib>
#include <Eigen/Core>
#include "Functions.h"
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class SuperGaussianModel : public Model
{
    public:
    
        SuperGaussianModel(const RefArrayXd covariates);
        ~SuperGaussianModel();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);

    protected:

    private:
    
}; // END class SuperGaussianModel


#endif
