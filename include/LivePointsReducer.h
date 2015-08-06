// Abstract class for nested sampling inference
// Enrico Corsaro @ IvS - 18 October 2013
// e-mail: emncorsaro@gmail.com
// Header file "LivePointsReducer.h"
// Implementation contained in "LivePointsReducer.cpp"

#ifndef LIVEPOINTSREDUCER_H
#define LIVEPOINTSREDUCER_H


#include <iostream>
#include <iomanip>
#include <cfloat>
#include <ctime>
#include <cmath>
#include <vector>
#include <cassert>
#include <limits>
#include <random>
#include <algorithm>
#include <Eigen/Dense>
#include "Functions.h"
#include "NestedSampler.h"


using namespace std;

class NestedSampler;

class LivePointsReducer
{
    public:

        LivePointsReducer(NestedSampler &nestedSampler);
        ~LivePointsReducer();
       
        vector<int> findIndicesOfLivePointsToRemove(mt19937 engine);
        int getNlivePointsToRemove();

        virtual int updateNlivePoints() = 0;
        

    protected:

        int NlivePointsAtCurrentIteration;
        int updatedNlivePoints;
        NestedSampler &nestedSampler;


    private:
        
        int NlivePointsToRemove;

};

#endif
