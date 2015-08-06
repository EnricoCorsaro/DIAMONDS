#include "PowerlawReducer.h"


// PowerlawReducer::PowerlawReducer()
//
// PURPOSE:
//      Derived class constructor. 
//
// INPUT:
//      nestedSampler:      a NestedSampler class object used as the container of
//                          information to use when reducing the number of live points.
//      exponent:           a double specifying the rate of the reduction process. For this 
//                          specific case, this number either enhances or smoothes the effect
//                          of the power law reduction. It is a number >= 0.
//                          Default is 1 meaning that a standard linear reduction occurs.
//                          If > 1 the reduction is suoer-linear, hence faster, and if
//                          < 1 the reduction is sub-linear, hence slower.
//      terminationFactor:  the fraction of remainder evidence to gained evidence used to terminate 
//                          the nested iteration loop. This value is also used as a tolerance on the final
//                          evidence to update the number of live points in the nesting process.
//

PowerlawReducer::PowerlawReducer(NestedSampler &nestedSampler, const double tolerance, 
                                       const double exponent, const double terminationFactor)
: LivePointsReducer(nestedSampler),
  tolerance(tolerance),
  exponent(exponent),
  terminationFactor(terminationFactor)
{
    assert(exponent >= 0.0);
    assert(tolerance >= 1.0);
}











// PowerlawReducer::~PowerlawReducer()
//
// PURPOSE:
//      Abstract base class destructor. 
//

PowerlawReducer::~PowerlawReducer()
{
}











// PowerlawReducer::updateNlivePoints()
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

int PowerlawReducer::updateNlivePoints()
{
    // Retrive the ratio of live to cumulated evidence for the current iteration 

    double ratioOfRemainderToCurrentEvidence = nestedSampler.getRatioOfRemainderToCurrentEvidence();        // Initial prior mass = 1


    // Evaluate the new number of live points to be used in the next iteration of the nesting process
   
    NlivePointsAtCurrentIteration = nestedSampler.getNlivePoints();
    double ratio = ratioOfRemainderToCurrentEvidence/terminationFactor;
    int NlivePointsToRemove = pow((tolerance/ratio), exponent);

    updatedNlivePoints = NlivePointsAtCurrentIteration - NlivePointsToRemove;


    // If new number of live points is lower than minNlivePoints, do not accept the new number and stick to the
    // previous one.

    if (updatedNlivePoints < nestedSampler.getMinNlivePoints())
    {
        updatedNlivePoints = NlivePointsAtCurrentIteration;
    }

    return updatedNlivePoints;
}
