
#include "RandomVariate.h"




RandomVariate::RandomVariate(int Ndim)
: Ndim(Ndim), minimum(Ndim), maximum(Ndim)
{
    minimum.setZero();
    maximum.setOnes();
}




RandomVariate::~RandomVariate()
{
    // pass
}







// RandomVariate::getNdim()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

int RandomVariate::getNdim()
{
    return Ndim;
}









// RandomVariate::setBoundaries()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:
// 
// REMARKS: FIXME: Better use the set() method of arrays?

void RandomVariate::setBoundaries(const RefArrayXd min, const RefArrayXd max)
{
    assert (min.size() == minimum.size());
    assert (max.size() == maximum.size());
    
    minimum = min;
    maximum = max;
}






// RandomVariate::setBoundaries()
//
// PURPOSE:
//      Overloaded setBoundaries() function, for convenience in the case Ndim = 1.
//
// INPUT:
//
// OUTPUT:

void RandomVariate::setBoundaries(double min, double max)
{
    assert (Ndim == 1);
    
    minimum(0) = min;
    maximum(0) = max;
}