
#include "RandomVariate.h"




RandomVariate::RandomVariate()
: minimum(0.0), maximum(1.0)
{
    // pass
}




RandomVariate::~RandomVariate()
{
    // pass
}




void RandomVariate::setBoundaries(double min, double max)
{
    minimum = min;
    maximum = max;
}