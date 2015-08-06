// Abstract class for nested sampling inference
// Enrico Corsaro @ IvS - 28 October 2013
// e-mail: emncorsaro@gmail.com
// Header file "PowerlawReducer.h"
// Implementation contained in "PowerlawReducer.cpp"

#ifndef POWERLAWREDUCER_H
#define POWERLAWREDUCER_H

#include "LivePointsReducer.h"

using namespace std;

class PowerlawReducer : public LivePointsReducer
{

    public:

        PowerlawReducer(NestedSampler &nestedSampler, const double tolerance = 100.0, 
                           const double exponent = 1, const double terminationFactor = 0.05);
        ~PowerlawReducer();
        
        virtual int updateNlivePoints();


    protected:

        double tolerance;
        double exponent;
        double terminationFactor;


    private:

};

#endif
