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
    
        PolynomialModel(const RefArrayXd covariates, const RefArrayXd covariatesUncertainties, const int Ndegrees, const double covariatesOffset);
        ~PolynomialModel();
        int getNdegrees();
        double getCovariatesOffset();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters);
        virtual void computeVariance(RefArrayXd modelVariance, const RefArrayXd modelParameters);

    protected:

        int Ndegrees;
        double covariatesOffset;
        ArrayXd covariatesUncertainties;

    private:
    
}; // END class PolynomialModel


#endif
