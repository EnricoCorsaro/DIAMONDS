// Derived class for building a zero model.
// Such a model always returns zero for any value of the fitparameters.


#ifndef ZEROMODEL_H
#define ZEROMODEL_H

#include <iostream>
#include "Model.h"


using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class ZeroModel : public Model
{
    public:
    
        ZeroModel(const RefArrayXd covariates);
        ~ZeroModel();

        virtual void predict(RefArrayXd predictions, const RefArrayXd modelParameters) override;
       
    protected:


    private:
    
};


#endif
