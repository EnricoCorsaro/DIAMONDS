#include "HyperQuadricSampler.h"

// HyperQuadricSampler::HyperQuadricSampler()
//
// PURPOSE: 
//      Base class constructor.
//
// INPUT:
//      prior: prior object containing the prior distributions for each parameter
//      likelihood: likelihood object containing the likelihood and the model to be used
//      metric: metric object to define the distance between two points
//      Nobjects: the initial number of objects coming from the main nested sampling loop
//

HyperQuadricSampler::HyperQuadricSampler(Prior &prior, Likelihood &likelihood, Metric &metric, const int Nobjects)
: Nobjects(Nobjects),
  prior(prior),
  likelihood(likelihood),
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
