// Abstract base class for building inference models
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 15 February 2013
// e-mail: emncorsaro@gmail.com
// Header file "Model.h"
// Implementations contained in "Model.cpp"


#ifndef MODEL_H
#define MODEL_H

#include <cstdlib>
#include <Eigen/Core>
#include "Functions.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class Model
{
    public:
    
        Model(const RefArrayXd covariates);
        ~Model();
        ArrayXd getCovariates();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters) = 0;
        int getNparameters();


    protected:
        
        ArrayXd covariates;
        int Nparameters;


    private:
    
}; // END class Model


#endif
