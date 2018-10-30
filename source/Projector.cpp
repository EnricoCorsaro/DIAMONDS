#include "Projector.h"

// Projector:Projector()
//
// PURPOSE: 
//      Class constructor
//
// INPUT:
//      printNdimensions: a boolean specifying whether one wants the reduced number
//      of dimensions to be explicitly printed on the screen.
//      Sets projected data member to 0 by default.
//

Projector::Projector(bool printNdimensions)
: printNdimensions(printNdimensions),
  reducedNdimensions(0)
{

}



// Projector::getReducedNdimensions()
//
// PURPOSE:
//      Get protected data member reducedNdimensions.
//
// OUTPUT:
//      An integer containing the final number of
//      nested loop iterations.
//

unsigned int Projector::getReducedNdimensions()
{
    return reducedNdimensions;
}
