// Abstract class for nested sampling inference
// Enrico Corsaro @ IvS - 28 October 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "ExponentialReducer.h"
// Implementation contained in "ExponentialReducer.cpp"

#ifndef EXPONENTIALREDUCER_H
#define EXPONENTIALREDUCER_H

#include "LivePointsReducer.h"

using namespace std;

class ExponentialReducer : public LivePointsReducer
{

    public:

        ExponentialReducer(NestedSampler &nestedSampler, const double tolerance = 10.0, 
                           const double enhancingFactor = 1, const double terminationFactor = 0.1);
        ~ExponentialReducer();
        
        virtual int updateNobjects();


    protected:

        double tolerance;
        double enhancingFactor;
        double terminationFactor;


    private:

};

#endif
