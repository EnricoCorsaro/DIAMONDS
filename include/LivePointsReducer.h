// Abstract class for nested sampling inference
// Enrico Corsaro @ IvS - 18 October 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
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
        int getNobjectsToRemove();

        virtual int updateNobjects() = 0;
        

    protected:

        int NobjectsAtCurrentIteration;
        int updatedNobjects;
        NestedSampler &nestedSampler;


    private:
        
        int NobjectsToRemove;

};

#endif
