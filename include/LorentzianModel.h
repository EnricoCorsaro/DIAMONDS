// Derived class for building mono lorentzian profile models
// Created by Enrico Corsaro @ IvS - 19 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "LorentzianModel.h"
// Implementations contained in "LorentzianModel.cpp"


#ifndef LORENTZIANMODEL_H
#define LORENTZIANMODEL_H

#include <iostream>
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;

class LorentzianModel : public Model
{
    public:
    
        LorentzianModel(const RefArrayXd covariates);
        ~LorentzianModel();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);
        int getNparameters();
       

    protected:


    private:

        int Nparameters;
    
}; // END class LorentzianModel


#endif
