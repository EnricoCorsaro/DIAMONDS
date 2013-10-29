// Abstract class for nested sampling inference
// Enrico Corsaro @ IvS - 18 October 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "LogarithmicReducer.h"
// Implementation contained in "LogarithmicReducer.cpp"

#ifndef LOGARITHMICREDUCER_H
#define LOGARITHMICREDUCER_H

#include "LivePointsReducer.h"

using namespace std;

class LogarithmicReducer : public LivePointsReducer
{

    public:

        LogarithmicReducer(NestedSampler &nestedSampler, const double reductionRate);
        ~LogarithmicReducer();
        
        virtual int updateNobjects();


    protected:

        double reductionRate;
        double informationGain;


    private:

};

#endif
