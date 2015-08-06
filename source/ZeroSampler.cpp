#include "ZeroSampler.h"

// ZeroSampler::ZeroSampler()
//
// PURPOSE: 
//      Constructor. Sets initial information, logEvidence and type 
//      of prior and likelihood distributions to be used. 
//
// INPUT:
//      printOnTheScreen:       Boolean value specifying whether the results are to 
//                              be printed on the screen or not.
//      initialNlivePoints:        Initial number of live points to start the nesting process
//      minNlivePoints:            Minimum number of live points allowed in the nesting process
//      ptrPriors:              Vector of pointers to Prior class objects
//      likelihood:             Likelihood class object used for likelihood sampling.
//      metric:                 Metric class object to contain the metric used in the problem.
//      clusterer:              Clusterer class object specifying the type of clustering algorithm to be used.
//

ZeroSampler::ZeroSampler(const bool printOnTheScreen, const int initialNlivePoints, const int minNlivePoints, vector<Prior*> ptrPriors, 
                    Likelihood &likelihood, Metric &metric, Clusterer &clusterer)
: NestedSampler(printOnTheScreen, initialNlivePoints, minNlivePoints, ptrPriors, likelihood, metric, clusterer)
{
}










// ZeroSampler::~ZeroSampler()
//
// PURPOSE: 
//      Base class destructor.
//

ZeroSampler::~ZeroSampler()
{

}











// ZeroSampler::drawWithConstraint()
//
// PURPOSE:
//      Empty function not to be used.
//
// INPUT:
//      totalSample:                Eigen Array matrix of size (Ndimensions, NlivePoints)
//                                  containing the total sample of active points at a given nesting iteration
//      Nclusters:                  Optimal number of clusters found by clustering algorithm
//      clusterIndices:             Indices of clusters for each point of the sample
//      clusterSizes:               A vector of integers containing the number of points belonging to each cluster
//      drawnPoint:                 Eigen Array matrix of size (Ndimensions,Ndraws) to contain the
//                                  coordinates of the drawn point to be used for the next nesting loop. 
//                                  When used for the first time, the array contains the coordinates of the 
//                                  worst nested object to be updated.
//      maxNdrawAttempts:           Maximum number of attempts allowed when drawing from a single ellipsoid.
//      
//
// OUTPUT:
//      A boolean value that is true if a new point in the sampling process is found and false otherwise.
//

bool ZeroSampler::drawWithConstraint(const RefArrayXXd totalSample, const unsigned int Nclusters, const vector<int> &clusterIndices,
                                               const vector<int> &clusterSizes, RefArrayXd drawnPoint, 
                                               double &logLikelihoodOfDrawnPoint, const int maxNdrawAttempts)
{    
    return false;
}
