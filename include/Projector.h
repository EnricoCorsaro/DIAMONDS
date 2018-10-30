// Abstract base class for implementing feature projection algorithms aimed
// at reducing the dimensionality of a given problem (e.g. inference problem).
// Created by Enrico Corsaro @ OACT - October 2018
// e-mail: emncorsaro@gmail.com
// Header file "Projector.h"
// Implementations contained in "Projector.cpp"


#ifndef PROJECTOR_H
#define PROJECTOR_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "Functions.h"

using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class Projector
{
    public:
    
        Projector(bool printNdimensions = false);
        ~Projector(){};
    
        virtual ArrayXXd projection(RefArrayXXd sample) = 0;
        unsigned int getReducedNdimensions();


    protected:
       
        bool printNdimensions;
        int Npoints;
        int Ndimensions;
        int reducedNdimensions;
    

    private:

};


#endif
