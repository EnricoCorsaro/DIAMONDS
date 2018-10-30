#include "Clusterer.h"

// Clusterer:Clusterer()
//
// PURPOSE:
//      Abstract base class constructor. 
//
// INPUT:
//      metric : the metric class to define a metric for the computation of distances between pair of points.
//      featureProjector: the class to evaluate a feature projection of the dataset to cluster, in order to
//                        reduce its dimensionality.
//      featureProjectionActivated: a boolean variable to decide whether or not performing the feature projection.
//

Clusterer::Clusterer(Metric &metric, Projector &featureProjector, bool featureProjectionActivated)
: metric(metric),
  featureProjector(featureProjector),
  featureProjectionActivated(featureProjectionActivated)
{

}










// Clusterer::getReducedNdimensions()
//
// PURPOSE:
//      Get protected data member reducedNdimensions of Projector class.
//
// OUTPUT:
//      An integer containing the final reduced dimensionality.
//

unsigned int Clusterer::getReducedNdimensions()
{
    return featureProjector.getReducedNdimensions();
}
