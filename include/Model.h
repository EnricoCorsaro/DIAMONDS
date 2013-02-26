// Abstract base class for building inference models
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 15 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "Model.h"
// Implementations contained in "Model.cpp"


#ifndef MODEL_H
#define MODEL_H

#include <cstdlib>
#include <Eigen/Core>
#include <MathExtra.h>


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class Model
{
    public:
    
        Model(const RefArrayXd covariates);
        ~Model();
        ArrayXd getCovariates();
        
        virtual void predict(RefArrayXd predictions, const RefArrayXd nestedSampleOfParameters) = 0;
        

    protected:
        
        ArrayXd covariates;


    private:
    
}; // END class Model


#endif
