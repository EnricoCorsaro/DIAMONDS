// Derived class for building a multi-dimensional Gaussian model
// Created by Enrico Corsaro @ OACT - December 2017
// e-mail: emncorsaro@gmail.com
// Header file "GaussianModel.h"
// Implementations contained in "GaussianModel.cpp"


#ifndef GAUSSIANMODEL_H
#define GAUSSIANMODEL_H

#include <cstdlib>
#include <Eigen/Core>
#include "Functions.h"
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class GaussianModel : public Model
{
    public:
    
        GaussianModel(const RefArrayXd covariates, int Nobservables);
        ~GaussianModel();
        ArrayXd getCovariates();
        int getNobservables();
        int getNpoints();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);

    protected:

        int Nobservables;
        int Npoints;

    private:
    
}; // END class GaussianModel


#endif
