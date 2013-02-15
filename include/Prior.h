
#ifndef PRIOR_H
#define PRIOR_H

#include <random>
#include <Eigen/Core>


using namespace std;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class Prior
{

    public:

        Prior(int Ndim);
        ~Prior();
        
        virtual void draw(RefArrayXXd parameterSample, const double Nobjects) = 0;
        virtual void drawWithConstraint(RefArrayXd parameter, Likelihood &likelihood) = 0;

    protected:

    private:
    
        int Ndim;

};

#endif
