
#ifndef UNIFORMPRIOR_H
#define UNIFORMPRIOR_H

#include <random>
#include <Eigen/Core>
#include "Prior.h"


using namespace std;
using Eigen:ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class UniformPrior : public Prior
{

    public:

        UniformPrior(const RefArrayXd min, const RefArrayXd max);
        ~UniformPrior();

        ArrayXd getMinimum();
        ArrayXd getMaximum();

        void draw(RefArrayXXd parameterSample, const double Nobjects);
        void drawWithConstraint(RefArrayXd parameter, Likelihood &likelihood);

    private:

        uniform_real_distribution<> uniform;
        mt19937 engine;
        ArrayXd min;
        ArrayXd max;

}; 

#endif
