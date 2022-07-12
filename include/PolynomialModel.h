// Derived class for building a polynomial model
// Created by Enrico Corsaro @ OACT - December 2017
// e-mail: emncorsaro@gmail.com
// Header file "PolynomialModel.h"
// Implementations contained in "PolynomialModel.cpp"


#ifndef POLYNOMIALMODEL_H
#define POLYNOMIALMODEL_H

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
    
        PolynomialModel(const RefArrayXd covariates, const int Ndegrees, const double covariatesOffset);
        ~PolynomialModel();
        int getNdegrees();
        double getCovariatesOffset();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);

    protected:

        int Ndegrees;
        double covariatesOffset;

    private:
    
}; // END class PolynomialModel


#endif
