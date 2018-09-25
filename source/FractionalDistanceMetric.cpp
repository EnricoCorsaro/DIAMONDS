#include "FractionalDistanceMetric.h"

// FractionalDistanceMetric::FractionalDistanceMetric()
//
// PURPOSE:
//      Class constructor.
//
// INPUT:
//      fraction:           The fraction (usually in the interval [0,1]) to compute the
//                          corresponding norm. The norm has the expression
//                          ( sum_i^d (x_i - y_i)^f )^(1/f), where f is the fraction term,
//                          and d is the number of dimensions of the dataset.
//                          See Aggarwal et al. "On the Surprising Behavior of Distance Metrics
//                          in High Dimensional Space" for more details.

FractionalDistanceMetric::FractionalDistanceMetric(const double fraction)
: fraction(fraction)
{
    assert(fraction > 0);
}





double FractionalDistanceMetric::distance(RefArrayXd point1, RefArrayXd point2)
{
    return pow((point1-point2).abs().pow(fraction).sum(),1.0/fraction);
}
