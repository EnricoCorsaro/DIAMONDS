// Abstract class for nested sampling inference
// Enrico Corsaro @ IvS - 18 October 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "FerozReducer.h"
// Implementation contained in "FerozReducer.cpp"

#ifndef FEROZREDUCER_H
#define FEROZREDUCER_H


#include "LivePointsReducer.h"


using namespace std;

class FerozReducer : public LivePointsReducer
{

    public:

        FerozReducer(NestedSampler &nestedSampler, const double toleranceOnEvidence);
        ~FerozReducer();
        
        virtual int updateNobjects();


    protected:

        double logMaxEvidenceContribution;          // The logarithm of the maximum evidence contribution at the previous iteration of the nesting process
        double toleranceOnEvidence;

    private:

};

#endif
