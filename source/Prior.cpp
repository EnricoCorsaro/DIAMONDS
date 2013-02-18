#include "Prior.h"

// Prior::Prior()
//
// PURPOSE:
//      Abstract base class constructor. Sets the dimensions
//      of the parameter space.
//
// INPUT:
//      Ndim : integer containing the number of dimensions, i.e.
//      the number of free parameters of the problem.
//
// OUTPUT:

Prior::Prior(const int Ndimensions, const int Nobjects)
: Ndimensions(Ndimensions), 
  Nobjects(Nobjects)
{

}







// Prior::~Prior()
//
// PURPOSE:
//      Abstract base class destructor. 
//
// INPUT:
/
// OUTPUT:

Prior::~Prior()
{

}







// Prior::getNdimensions()
//
// PURPOSE: 
//      Get function to obtain the protected data member Ndimensions.
//
// INPUT:
//
// OUTPUT:

int Prior::getNdimensions()
{
    return Ndimensions;
}







// Prior::getNobjects()
//
// PURPOSE: 
//      Get function to obtained the protected data member Nobjects.
//
// INPUT:
//
// OUTPUT:

int Prior::getNobjects()
{
    return Nobjects;
}
