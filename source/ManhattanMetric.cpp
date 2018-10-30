#include "ManhattanMetric.h"


double ManhattanMetric::distance(RefArrayXd point1, RefArrayXd point2)
{
    return (point1-point2).abs().sum();
}
