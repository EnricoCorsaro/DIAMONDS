// Derived class for building a multilinear model
// Created by Enrico Corsaro @ OACT - December 2017
// e-mail: emncorsaro@gmail.com
// Header file "MultiLinearModel.h"
// Implementations contained in "MultiLinearModel.cpp"


#ifndef MULTILINEARMODEL_H
#define MULTILINEARMODEL_H

#include <cstdlib>
#include <Eigen/Core>
#include "Functions.h"
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class MultiLinearModel : public Model
{
    public:
    
        MultiLinearModel(const RefArrayXd covariates, const RefArrayXd covariatesUncertainties, int Nobservables);
        ~MultiLinearModel();
        int getNobservables();
        ArrayXd getCovariatesUncertainties();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);
        virtual void computeVariance(RefArrayXd modelVariance, const RefArrayXd modelParameters);

    protected:

        int Nobservables;
        ArrayXd covariatesUncertainties;


    private:
    
}; // END class MultiLinearModel


#endif
