#include "HyperQuadricSampler.h"

// HyperQuadricSampler::HyperQuadricSampler()
//
// PURPOSE: 
//      Base class constructor.
//

HyperQuadricSampler::HyperQuadricSampler(Prior &prior, Metric &metric, const int Nobjects)
: Nobjects(Nobjects),
  prior(prior),
  metric(metric)
{

}











// HyperQuadricSampler::~HyperQuadricSampler()
//
// PURPOSE: 
//      Base class destructor.
//

HyperQuadricSampler::~HyperQuadricSampler()
{

}
