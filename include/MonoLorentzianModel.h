// Derived class for building mono lorentzian profile models
// Created by Enrico Corsaro @ IvS - 19 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "MonoLorentzianModel.h"
// Implementations contained in "MonoLorentzianModel.cpp"


#ifndef MONOLORENTZIANMODEL_H
#define MONOLORENTZIANMODEL_H

#include <iostream>
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class MonoLorentzianModel
{
    public:
    
        MonoLorentzianModel(const RefArrayXd covariates);
        ~MonoLorentzianModel();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);
        int getParametersNumber();
        
    protected:


    private:

        int parametersNumber;
    
}; // END class MonoLorentzianModel


#endif
