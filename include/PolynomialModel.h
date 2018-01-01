// Derived class for building a linear model
// Created by Enrico Corsaro @ OACT - December 2017
// e-mail: emncorsaro@gmail.com
// Header file "PolynomialModel.h"
// Implementations contained in "PolynomialModel.cpp"


#ifndef LINEARMODEL_H
#define LINEARMODEL_H

#include <cstdlib>
#include <Eigen/Core>
#include "Functions.h"
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class PolynomialModel : public Model
{
    public:
    
        PolynomialModel(const RefArrayXd covariates, int Ndegrees);
        ~PolynomialModel();
        ArrayXd getCovariates();
        int getNdegrees();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);

    protected:

        int Ndegrees;

    private:
    
}; // END class PolynomialModel


#endif
