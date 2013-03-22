
#ifndef METRIC_H
#define METRIC_H


#include <Eigen/Core>

using namespace std;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class Metric
{
    public:
    
        Metric(){};
        ~Metric(){};

        virtual double distance(RefArrayXd point1, RefArrayXd point2) = 0;

    protected:
    
    private:
};


#endif