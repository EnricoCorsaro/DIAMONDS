// Class for the computation of the
// theoretical models.
// Enrico Corsaro @ IvS - 14 February 2013
// e-mail: enrico.corsaro@ster.kuleuven.be
// Header file "Model.h"
// Implementation contained in "Model.cpp"

#ifndef MODEL_H
#define MODEL_H

#include <Eigen/Core>
#include "MathExtra.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;

class Model
{

    public:
        
        Model(const int modelID);
        ~Model();
        void setModel(const RefArrayXd x_obs, const RefArrayXd setOfParameters);
        int getParametersNumber();
        ArrayXd getYTheor();

    private:

        int modelIdentifier;
        int parametersNumber;
        ArrayXd x;
        ArrayXd y_theor;
        
        void buildModel1();
        void buildModel2();
        void buildModel3();

}; // END class Model

#endif
