#include "FerozReducer.h"


// FerozReducer::FerozReducer()
//
// PURPOSE:
//      Derived class constructor. 
//
// INPUT:
//      nestedSampler:  a NestedSampler class object used as the container of
//                      information to use when reducing the number of live points.
//

FerozReducer::FerozReducer(NestedSampler &nestedSampler, const double toleranceOnEvidence)
: LivePointsReducer(nestedSampler),
  toleranceOnEvidence(toleranceOnEvidence)
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













// FerozReducer::updateNobjects()
//
// PURPOSE:
//      Updates the number of live points for the next iteration of the nesting process.
//      For this case, the equation adopted is that used by
//      Feroz F. et al. (2009), MNRAS, 398, 1601.
//
// OUTPUT:
//      An integer specifying the final number of live points to be adopted.
//

int FerozReducer::updateNobjects()
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
   
    NobjectsAtCurrentIteration = nestedSampler.getNobjects();
    double numerator = exp(Functions::logExpDifference(logMaxEvidenceContribution,logMaxEvidenceContributionNew));
    double denominator = exp(logMaxEvidenceContributionNew);
    updatedNobjects = NobjectsAtCurrentIteration - static_cast<int>(nestedSampler.getMinNobjects() * numerator / (denominator * toleranceOnEvidence));
  

    // If new number of live points is lower than zero, do not accept the new number and stick to the
    // previous one.

    if (updatedNobjects < 0)
    {
        updatedNobjects = NobjectsAtCurrentIteration;
    }


    // Finally update max evidence contribution with newest value
    
    logMaxEvidenceContribution = logMaxEvidenceContributionNew;

    return updatedNobjects;
}

