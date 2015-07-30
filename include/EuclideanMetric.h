// Derived class for creating an euclidean metric objects
// used within the clustering class.
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 11 April 2013
// e-mail: emncorsaro@gmail.com
// Header file "EuclideanMetric.h"
// Implementation contained in "EuclideanMetric.cpp"

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
