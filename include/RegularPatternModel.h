// Derived class for building mono lorentzian profile models
// Created by Enrico Corsaro @ IvS - 13 June 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "RegularPatternModel.h"
// Implementations contained in "RegularPatternModel.cpp"


#ifndef REGULARPATTERNMODEL_H
#define REGULARPATTERNMODEL_H

#include <iostream>
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class RegularPatternModel : public Model
{
    public:
    
        RegularPatternModel(const RefArrayXd covariates, const int Norders);
        ~RegularPatternModel();

        virtual void predict(RefArrayXd predictions, const RefArrayXd nestedSampleOfParameters);
        int getNparameters();
       

    protected:


    private:

        int Norders;
        int Nparameters;
   

}; // END class RegularPatternModel


#endif
