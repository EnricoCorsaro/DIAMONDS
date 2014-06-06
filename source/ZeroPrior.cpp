#include "ZeroPrior.h"



// ZeroPrior::ZeroPrior()
//
// PURPOSE: 
//      Derived class constructor.
//
// INPUT:
//

ZeroPrior::ZeroPrior(const int Ndimensions)
: Prior(Ndimensions)
{
}









// ZeroPrior:~ZeroPrior()
//
// PURPOSE: 
//      Derived class destructor.
//

ZeroPrior::~ZeroPrior()
{

}










// ZeroPrior::logDensity()
//
// PURPOSE:
//      Returns a -Infinity logDensity for any input coordinate.
//
// INPUT: 
//      x:                   Point in which the log(pdf) should be evaluated.
//      includeConstantTerm: If true : compute the exact log(density), 
//                           If false: ignore the constant terms (with factors of pi, 2, etc.)
//
// OUTPUT:
//      Natural logarithm of the probability density evaluation in x.
//

double ZeroPrior::logDensity(RefArrayXd const x, const bool includeConstantTerm)
{
    double logDens = minusInfinity;

    return logDens;
}













// ZeroPrior::drawnPointIsAccepted()
//
// PURPOSE:
//      Return false for any input drawnPoint coordinates
//
// INPUT: 
//      drawnPoint:     an Eigen array containing the coordinates 
//                      of the new drawn point to be verified
//
// OUTPUT:
//      Returns true if the new point is accepted, false if not.
//

bool ZeroPrior::drawnPointIsAccepted(RefArrayXd const drawnPoint)
{
    return false;
}
