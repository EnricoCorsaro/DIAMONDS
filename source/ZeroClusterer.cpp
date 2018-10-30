#include "ZeroClusterer.h"


// ZeroClusterer::ZeroClusterer()
//
// PURPOSE: 
//      Constructor.
//
// INPUT: 
//      metric: this class is used to compute the squared distance between two points
// 

ZeroClusterer::ZeroClusterer(Metric &metric, Projector &featureProjector, bool featureProjectionActivated)
: Clusterer(metric, featureProjector, featureProjectionActivated) 
{
}









// ZeroClusterer::~ZeroClusterer()
//
// PURPOSE: 
//      Destructor.
//

ZeroClusterer::~ZeroClusterer()
{

}









// ZeroClusterer::cluster()
//
// PURPOSE: 
//      Empty function, not used when building a ZeroClusterer object
//
// INPUT:
//      sample(Ndimensions, Npoints): sample of N-dimensional points
//      optimalClusterIndices(Npoints): for each point the index of the cluster it belongs to. This index
//                                      runs from 0 to Nclusters-1.
//      optimalClusterSizes(Nclusters): for each of the clusters, this vector contains the number of points
// 
// OUTPUT:
//      An integer equal to zero.
//

int ZeroClusterer::cluster(RefArrayXXd sample, vector<int> &optimalClusterIndices, vector<int> &optimalClusterSizes)
{
    int zeroValue = 0;
    return zeroValue;
}
