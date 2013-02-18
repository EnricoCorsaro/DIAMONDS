
#ifndef UNIFORMPRIOR_H
#define UNIFORMPRIOR_H

#include <cassert>
#include <iostream>
#include <Eigen/Core>
#include "Prior.h"


using namespace std;
using Eigen:ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class UniformPrior : public Prior
{

    public:

        UniformPrior(const RefArrayXd minimum, const RefArrayXd maximum), const int Nobjects;
        ~UniformPrior();

        ArrayXd getMinimum();
        ArrayXd getMaximum();
        double getUniformFactor();

        virtual void draw(RefArrayXXd nestedParameters);
        virtual void drawWithConstraint(RefArrayXd nestedParameters, Likelihood &likelihood);

    private:

        uniform_real_distribution<> uniform;
        mt19937 engine;
        ArrayXd minimum;
        ArrayXd maximum;
        double uniformFactor;

}; 

#endif
