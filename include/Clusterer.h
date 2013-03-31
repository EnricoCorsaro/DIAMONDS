
#ifndef CLUSTERER_H
#define CLUSTERER_H

#include <vector>
#include <Eigen/Core>
#include "Metric.h"


using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;
typedef Eigen::Ref<Eigen::ArrayXi> RefArrayXi;


class Clusterer
{
    public:
    
        Clusterer(Metric &metric);
        ~Clusterer(){};    
        virtual int cluster(RefArrayXXd sample, RefArrayXi clusterIndices) = 0;


    protected:
        
        Metric &metric;
    

    private:

};




#endif
