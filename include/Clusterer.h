
#ifndef CLUSTERER_H
#define CLUSTERER_H

#include <vector>
#include <Eigen/Core>
#include "Metric.h"


using namespace std;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class Clusterer
{
    public:
    
        Clusterer(Metric &metric);
        ~Clusterer(){};    
        virtual int cluster(RefArrayXXd sample, vector<int> &clusterIndices) = 0;

        Metric &metric;

    protected:
    
    private:

};




#endif