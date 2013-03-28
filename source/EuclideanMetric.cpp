
#include "EuclideanMetric.h"


double EuclideanMetric::distance(RefArrayXd point1, RefArrayXd point2)
{
    return sqrt((point1-point2).square().sum());
}