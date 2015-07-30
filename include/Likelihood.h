// Abstract base class for likelihood computations.
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 15 February 2013
// e-mail: emncorsaro@gmail.com
// Header file "Likelihood.h"
// Implementations contained in "Likelihood.cpp"

#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <Eigen/Core>
#include "Functions.h"
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class Likelihood
{

    public:

        Likelihood(const RefArrayXd observations, Model &model);
        ~Likelihood();
        ArrayXd getObservations();

        virtual double logValue(RefArrayXd const modelParameters) = 0;


    protected:
        
        ArrayXd observations;
        Model &model;


    private:

        
};

#endif
