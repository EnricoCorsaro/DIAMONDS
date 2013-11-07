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

        ExponentialReducer(NestedSampler &nestedSampler, const double enhancingFactor);
        ~ExponentialReducer();
        
        virtual int updateNobjects();


    protected:

        double enhancingFactor;
        double informationGain;


    private:

};

#endif
