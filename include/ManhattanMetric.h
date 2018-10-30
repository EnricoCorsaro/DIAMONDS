// Derived class for creating an euclidean metric objects
// used within the clustering class.
// Created by Enrico Corsaro @ OACT - September 2018
// e-mail: emncorsaro@gmail.com
// Header file "ManhattanMetric.h"
// Implementation contained in "ManhattanMetric.cpp"

#ifndef MANHATTANMETRIC_H
#define MANHATTANMETRIC_H


#include "Metric.h"


using namespace std;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class ManhattanMetric : public Metric
{
    public:
    
        ManhattanMetric(){};
        ~ManhattanMetric(){};

        virtual double distance(RefArrayXd point1, RefArrayXd point2);

    protected:
    
    private:
    
};


#endif
