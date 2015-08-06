#include "LivePointsReducer.h"


// LivePointsReducer::LivePointsReducer()
//
// PURPOSE:
//      Abstract base class constructor. 
//
// INPUT:
//
//      nestedSampler:  a NestedSampler class object used as the container of
//                      information to use when reducing the number of live points.
//

LivePointsReducer::LivePointsReducer(NestedSampler &nestedSampler)
: nestedSampler(nestedSampler)
{
}











// LivePointsReducer::~LivePointsReducer()
//
// PURPOSE:
//      Abstract base class destructor. 
//

LivePointsReducer::~LivePointsReducer()
{
}












// LivePointsReducer::findIndicesOfLivePointsToRemove()
//
// PURPOSE:
//          Picks one live point randomly and save its index. 
//          Repeats the process up to the total number of live points to be removed.
//
// INPUT:
//          engine:     a Marsenne-Twister engine containing the clock seed for
//                      generating random numbers
//
// OUTPUT:
//          A vector to contain the indices of the live points to be removed.
//

vector<int> LivePointsReducer::findIndicesOfLivePointsToRemove(mt19937 engine)
{
    // Compute how many live points must be removed from the current sample.
    // In case no live points must be removed, process will skip initialization
    // of indices.
    
    NlivePointsToRemove = NlivePointsAtCurrentIteration - updatedNlivePoints;


    // Create an array of integers to contain the indices of the live points to remove

    vector<int> indicesOfLivePointsToRemove;

    if (NlivePointsToRemove > 0)
    {
        // At least one live point has to be removed, hence proceed with inizializing the
        // vector of indices

        for (int m = 0; m < NlivePointsToRemove; ++m)
        {
            // Create uniform random integer generator with current number of live points.
            // This is done because we remove one live point per time. Everytime we repeat the
            // process, we have to take into account that the total number of live points
            // was descreased by one.

            uniform_int_distribution<int> discreteUniform(0, NlivePointsAtCurrentIteration-1);
                

            // Select randomly one live point from the actual sample

            indicesOfLivePointsToRemove.push_back(discreteUniform(engine));
               

            // Reduce the current number of live points by one.
                
            --NlivePointsAtCurrentIteration;
        }
    }


    // If the conditional statement is skipped, the returned vector has zero elements.

    return indicesOfLivePointsToRemove;
}











// LivePointsReducer::getNlivePointsToRemove()
//
// PURPOSE:
//      Gets private data member NlivePointsToRemove
//

int LivePointsReducer::getNlivePointsToRemove()
{
    return NlivePointsToRemove;
}
