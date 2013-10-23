#include "LogarithmicReducer.h"


// LogarithmicReducer::LogarithmicReducer()
//
// PURPOSE:
//      Derived class constructor. 
//
// INPUT:
//      nestedSampler:  a NestedSampler class object used as the container of
//                      information to use when reducing the number of live points.
//

LogarithmicReducer::LogarithmicReducer(NestedSampler &nestedSampler)
: LivePointsReducer(nestedSampler)
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
