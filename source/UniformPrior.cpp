
#include "UniformPrior.h"



// UniformPrior:UniformPrior()
//
// PURPOSE: constructor
//
// INPUT:
//      min:
//      max:
// 
// OUTPUT:
//

UniformPrior::UniformPrior(const RefArrayXd min, const RefArrayXd max)
: Prior(min.size()), min(min), max(max)
{
    assert (min.size() == max.size());
}





// UniformPrior:~UniformPrior()
//
// PURPOSE: destructor
//

UniformPrior::~UniformPrior()
{

}


