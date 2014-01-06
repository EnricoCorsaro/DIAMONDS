#include "ExponentialReducer.h"


// ExponentialReducer::ExponentialReducer()
//
// PURPOSE:
//      Derived class constructor. 
//
// INPUT:
//      nestedSampler:      a NestedSampler class object used as the container of
//                          information to use when reducing the number of live points.
//      enhancingFactor:    a double specifying the rate of the reduction process. For this 
//                          specific case, this number either enhances or smoothes the effect
//                          of the exponential reduction. It is a number > 0.
//                          Default is 1 meaning that a standard exponential reduction occurs.
//      
//

ExponentialReducer::ExponentialReducer(NestedSampler &nestedSampler, const double tolerance, 
                                       const double enhancingFactor, const double terminationFactor)
: LivePointsReducer(nestedSampler),
  tolerance(tolerance),
  enhancingFactor(enhancingFactor),
  terminationFactor(terminationFactor)
{
    assert(enhancingFactor >= 0.0);
    assert(tolerance >= 1.0);
}











// ExponentialReducer::~ExponentialReducer()
//
// PURPOSE:
//      Abstract base class destructor. 
//

ExponentialReducer::~ExponentialReducer()
{
}











// ExponentialReducer::updateNobjects()
//
// PURPOSE:
//      Updates the number of live points for the upcoming iteration of the nesting process.
//      For this case, the exponential expression adopted allows to start reducing live points
//      after the tolerance on the ratio of the live to the cumulated evidence has been reached. 
//
// OUTPUT:
//      An integer specifying the final number of live points to be adopted.
//
// REMARK:
//      The returned value of live points is ensured to be not below the minimum allowed. 
//

int ExponentialReducer::updateNobjects()
{
    // Retrive the ratio of live to cumulated evidence for the current iteration 

    double ratioOfRemainderToCurrentEvidence = nestedSampler.getRatioOfRemainderToCurrentEvidence();        // Initial prior mass = 1


    // Evaluate the new number of live points to be used in the next iteration of the nesting process
   
    NobjectsAtCurrentIteration = nestedSampler.getNobjects();
    double ratio = ratioOfRemainderToCurrentEvidence/terminationFactor;
    double exponent = -1.0*enhancingFactor*(ratio - tolerance);
    int NobjectsToRemove = exp(exponent);
    updatedNobjects = NobjectsAtCurrentIteration - NobjectsToRemove;


    // If new number of live points is lower than minNobjects, do not accept the new number and stick to the
    // previous one.

    if (updatedNobjects < nestedSampler.getMinNobjects())
    {
        updatedNobjects = NobjectsAtCurrentIteration;
    }

    return updatedNobjects;
}
