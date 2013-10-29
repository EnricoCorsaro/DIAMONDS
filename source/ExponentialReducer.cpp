#include "ExponentialReducer.h"


// ExponentialReducer::ExponentialReducer()
//
// PURPOSE:
//      Derived class constructor. 
//
// INPUT:
//      nestedSampler:  a NestedSampler class object used as the container of
//                      information to use when reducing the number of live points.
//      reductionRate:  a double specifying the rate of the reduction process.
//

ExponentialReducer::ExponentialReducer(NestedSampler &nestedSampler, const double reductionRate)
: LivePointsReducer(nestedSampler),
  reductionRate(reductionRate)
{
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
//      Updates the number of live points for the next iteration of the nesting process.
//      For this case, the exponential expression adopted allows to reduce less points at the beginning
//      and speeds up while the number of live points approaches the minimum allowed. 
//
// OUTPUT:
//      An integer specifying the final number of live points to be adopted.
//

int ExponentialReducer::updateNobjects()
{
    if (nestedSampler.getNiterations() == 0)
    {
        // For the particular case of the first iteration initialize informationGain at the beginning.
   
        informationGain = 0.0;
    }


    // Evaluate new informationGain for the current iteration 

    double informationGainNew = nestedSampler.getInformationGain();


    // Evaluate the new number of live points to be used in the next iteration of the nesting process
   
    NobjectsAtCurrentIteration = nestedSampler.getNobjects();
    double numerator = exp(-1.0*(NobjectsAtCurrentIteration - nestedSampler.getMinNobjects()))
                       - exp(-1.0*(nestedSampler.getInitialNobjects() - nestedSampler.getMinNobjects())) 
                       * exp(informationGain - informationGainNew);
    double denominator = 1 - exp(-1.0*(nestedSampler.getInitialNobjects() - nestedSampler.getMinNobjects()));
    updatedNobjects = NobjectsAtCurrentIteration - static_cast<int>(nestedSampler.getMinNobjects() * (numerator * reductionRate) / denominator);


    // If new number of live points is lower than minNobjects, do not accept the new number 
    // and stick to the previous one.

    if (updatedNobjects < nestedSampler.getMinNobjects())
    {
        updatedNobjects = NobjectsAtCurrentIteration;
    }


    // Finally update information gain with newest value
    
    informationGain = informationGainNew;

    return updatedNobjects;
}
