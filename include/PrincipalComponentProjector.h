// Abstract base class for implementing feature projection algorithms aimed
// at reducing the dimensionality of a given problem (e.g. inference problem).
// Created by Enrico Corsaro @ OACT - October 2018
// e-mail: emncorsaro@gmail.com
// Header file "PrincipalComponentProjector.h"
// Implementations contained in "PrincipalComponentProjector.cpp"


#ifndef PRINCIPALCOMPONENTPROJECTOR_H
#define PRINCIPALCOMPONENTPROJECTOR_H

#include "Projector.h"


using namespace std;
using namespace Eigen;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;


class PrincipalComponentProjector : public Projector
{
    public:
    
        PrincipalComponentProjector(bool printNdimensions, double scalingFactor = 0.7);
        ~PrincipalComponentProjector(){};
    
        virtual ArrayXXd projection(RefArrayXXd sample);


    protected:
    
        double scalingFactor;
    

    private:

};


#endif
