#include "LogarithmicReducer.h"


// LogarithmicReducer::LogarithmicReducer()
//
// PURPOSE:
//      Derived class constructor. 
//
// INPUT:
//      nestedSampler:  a NestedSampler class object used as the container of
//                      information to use when reducing the number of live points.
//      reductionRate:  a double specifying the rate of the reduction process. This number
//                      has to be > 0. For the case = 1, no change in the reduction function
//                      is done. For accelerating the reduction process, the reduction
//                      rate must be increased.
//

LogarithmicReducer::LogarithmicReducer(NestedSampler &nestedSampler, const double reductionRate)
: LivePointsReducer(nestedSampler),
  reductionRate(reductionRate)
{
}











// LogarithmicReducer::~LogarithmicReducer()
//
// PURPOSE:
//      Abstract base class destructor. 
//

LogarithmicReducer::~LogarithmicReducer()
{
}











// LogarithmicReducer::updateNobjects()
//
// PURPOSE:
//      Updates the number of live points for the next iteration of the nesting process.
//      For this case, the logarithmic expression adopted allows to reduce more points at the beginning
//      and slows down while the number of live points approaches the minimum allowed. 
//
// OUTPUT:
//      An integer specifying the final number of live points to be adopted.
//

int LogarithmicReducer::updateNobjects()
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
    double numerator = log(nestedSampler.getInitialNobjects() - NobjectsAtCurrentIteration + 2) * exp(informationGain - informationGainNew);
    double denominator = log(nestedSampler.getInitialNobjects() - nestedSampler.getMinNobjects() + 2);
    updatedNobjects = NobjectsAtCurrentIteration - static_cast<int>(nestedSampler.getMinNobjects() * (numerator * reductionRate) / denominator);


    // Finally update information gain with newest value
    
    informationGain = informationGainNew;

    return updatedNobjects;
}
