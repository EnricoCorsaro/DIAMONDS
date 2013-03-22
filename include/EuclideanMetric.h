#ifndef EUCLIDEANMETRIC_H
#define EUCLIDEANMETRIC_H


#include <Eigen/Core>
#include "Metric.h"


using namespace std;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class EuclideanMetric : public Metric
{
    public:
    
        EuclideanMetric(){};
        ~EuclideanMetric(){};

        virtual double distance(RefArrayXd point1, RefArrayXd point2);

    protected:
    
    private:
    
};


#endif