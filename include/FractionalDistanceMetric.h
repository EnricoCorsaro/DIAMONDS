// Derived class for creating an euclidean metric objects
// used within the clustering class.
// Created by Enrico Corsaro @ OACT    - September 2018
// e-mail: emncorsaro@gmail.com
// Header file "FractionalDistanceMetric.h"
// Implementation contained in "FractionalDistanceMetric.cpp"

#ifndef FRACTIONALDISTANCEMETRIC_H
#define FRACTIONALDISTANCEMETRIC_H


#include "Metric.h"
#include <cassert>


using namespace std;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class FractionalDistanceMetric : public Metric
{
    public:
    
        FractionalDistanceMetric(const double fraction = 0.3);
        ~FractionalDistanceMetric(){};

        virtual double distance(RefArrayXd point1, RefArrayXd point2);

    protected:

        double fraction;
    
    private:
    
};


#endif
