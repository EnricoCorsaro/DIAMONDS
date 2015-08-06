#include "FerozReducer.h"


// FerozReducer::FerozReducer()
//
// PURPOSE:
//      Derived class constructor. 
//
// INPUT:
//      nestedSampler:  a NestedSampler class object used as the container of
//                      information to use when reducing the number of live points.
//      tolerance:      a double containing the tolerance on the evidence. The lower the factor
//                      the more live points are removed at each iteration.
//

FerozReducer::FerozReducer(NestedSampler &nestedSampler, const double tolerance)
: LivePointsReducer(nestedSampler),
  tolerance(tolerance)
{
}










// FerozReducer::~FerozReducer()
//
// PURPOSE:
//      Abstract base class destructor. 
//

FerozReducer::~FerozReducer()
{
}













// FerozReducer::updateNlivePoints()
//
// PURPOSE:
//      Updates the number of live points for the next iteration of the nesting process.
//      For this case, the equation adopted is that used by
//      Feroz F. et al. (2009), MNRAS, 398, 1601.
//
// OUTPUT:
//      An integer specifying the final number of live points to be adopted.
//
// REMARK:
//      The returned value of live points is ensured to be not below the minimum allowed. 
//

int FerozReducer::updateNlivePoints()
{
    if (nestedSampler.getNiterations() == 0)
    {
        // For the particular case of the first iteration initialize logMaxEvidenceContribution at the beginning.
        // Evaluate max evidence contribution for first iteration based on the logarithm of the 
        // maximum likelihood value of the initial set of live points
    
        logMaxEvidenceContribution = nestedSampler.getLogMaxLikelihoodOfLivePoints();        // Initial prior mass = 1
    }
    

    // Evaluate max evidence contribution for the current iteration 

    double logMaxEvidenceContributionNew = nestedSampler.getLogLikelihood().maxCoeff() + nestedSampler.getLogRemainingPriorMass();


    // Evaluate the new number of live points to be used in the next iteration of the nesting process
   
    NlivePointsAtCurrentIteration = nestedSampler.getNlivePoints();
    double numerator = exp(Functions::logExpDifference(logMaxEvidenceContribution,logMaxEvidenceContributionNew));
    double denominator = exp(logMaxEvidenceContributionNew);
    int NlivePointsToRemove = nestedSampler.getMinNlivePoints() * numerator / (denominator * tolerance);
    updatedNlivePoints = NlivePointsAtCurrentIteration - NlivePointsToRemove;
    

    // If new number of live points is lower than minNlivePoints, do not accept the new number and stick to the
    // previous one.

    if (updatedNlivePoints < nestedSampler.getMinNlivePoints())
    {
        updatedNlivePoints = NlivePointsAtCurrentIteration;
    }


    // Finally update max evidence contribution with newest value
    
    logMaxEvidenceContribution = logMaxEvidenceContributionNew;

    return updatedNlivePoints;
}

