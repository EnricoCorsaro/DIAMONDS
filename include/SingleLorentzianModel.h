// Derived class for building single lorentzian profile models
// Created by Enrico Corsaro @ IvS - 19 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "SingleLorentzianModel.h"
// Implementations contained in "SingleLorentzianModel.cpp"


#ifndef SINGLELORENTZIANMODEL_H
#define SINGLELORENTZIANMODEL_H

#include <iostream>
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class SingleLorentzianModel
{
    public:
    
        SingleLorentzianModel(const RefArrayXd covariates);
        ~SingleLorentzianModel();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);
        int getParametersNumber();
        
    protected:


    private:

        int parametersNumber;
    
}; // END class SingleLorentzianModel


#endif
