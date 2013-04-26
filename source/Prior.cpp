#include "Prior.h"

// Prior::Prior()
//
// PURPOSE:
//      Abstract base class constructor. Sets the dimensions
//      of the parameter space.
//
// INPUT:
//      Ndimensions : integer containing the number of dimensions, i.e.
//      the number of free parameters of the problem.
//      uniformFlag: a boolean variable specifying whether the prior
//      is uniform or not.
//

Prior::Prior(const int Ndimensions, const bool uniformFlag)
: Ndimensions(Ndimensions),
  uniformFlag(uniformFlag),
  engine(time(0))
{

}









// Prior::~Prior()
//
// PURPOSE:
//      Abstract base class destructor. 
//

Prior::~Prior()
{

}










// Prior::priorIsUniform()
//
// PURPOSE: 
//      States whether the prior is a uniform prior or not.
//
// OUTPUT:
//      A boolean variable specifying whether the prior is uniform (true)
//      or different than uniform (false).

bool Prior::priorIsUniform()
{
    return uniformFlag;
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


