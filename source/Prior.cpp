#include "Prior.h"

// Prior::Prior()
//
// PURPOSE:
//      Abstract base class constructor. Sets the dimensions
//      of the parameter space.
//
// INPUT:
//      Ndimensions : integer containing the number of dimensions, i.e.
//                    the number of free parameters of the problem.
//

Prior::Prior(const int Ndimensions)
: minusInfinity(numeric_limits<double>::lowest()),
  Ndimensions(Ndimensions)
{
	// Set the seed of the random generator using the clock

    clock_t clockticks = clock();
    engine.seed(clockticks);
}









// Prior::~Prior()
//
// PURPOSE:
//      Abstract base class destructor. 
//

Prior::~Prior()
{

}









// Prior::getNdimensions()
//
// PURPOSE: 
//      Get function to obtain the protected data member Ndimensions.
//
// OUTPUT:
//      An integer containing the number of dimensions of the parameter space.
//

int Prior::getNdimensions()
{
    return Ndimensions;
}


