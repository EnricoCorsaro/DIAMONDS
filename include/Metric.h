// Class for creating an abstract class to define metric distances
// used within the sampler class.
// Created by Enrico Corsaro & Joris De Ridder @ IvS - 11 April 2013
// e-mail: emncorsaro@gmail.com
// Header file "Metric.h"
// Implementation contained in "Metric.cpp"

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
