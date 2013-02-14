#include "Likelihood.h"

// Likelihood::Likelihood()
//
// PURPOSE:
//      Constructor. Reads input data file and chooses 
//      the theoretical model to build the likelihood.
//
// INPUT:
//
// OUTPUT:

Likelihood::Likelihood(const RefArrayXXd data, const int modelID)
: modelIdentifier(modelID)
{
    int Ncols = 0;

    Ncols = data.row(0).size();

    if (Ncols > 3)
    {
        cerr << "The number of columns of the input data file exceeds the format.\nQuitting program." << endl;
        exit(1);
    }
    else if (Ncols == 2)
    {
        covariates = data.col(0);
        observations = data.col(1);
        uncertainties = 0;
    }
    covariates = data.col(0);
    observations = data.col(1);
    uncertainties = data.col(2);

}






// Likelihood::setLikelihood()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

void Likelihood::setLikelihood(const RefArrayXd modelParameters)
{
    Model model(covariates, modelParameters, modelIdentifier);
    predictions = model.getPredictions();

    return;
}







// Likelihood::buildLogGaussianLikelihood()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

double Likelihood::buildLogGaussianLikelihood()
{
    double sumSigma = 0;

    sumSigma = MathExtra::sum(uncertainties);
    
    if (sigmSigma == 0.)
    {
        cerr << "Gaussian Likelihood cannot be computed. Missing uncertainty information." 
        << endl << "Quitting program." << endl;
        exit(1);
    }
    

    double likelihoodValue;
    likelihoodValue = MathExtra::logGaussLikelihood(observations, predictions, uncertainties);

    return likelihoodValue;
}







// Likelihood::buildMarginalizedLogGaussianLikelihood()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

double Likelihood::buildMarginalizedLogGaussianLikelihood()
{
    double likelihoodValue;

    return likelihoodValue;
}








// Likelihood::buildMedianLikelihood()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

double Likelihood::buildMedianLikelihood()
{
    double likelihoodValue;

    return likelihoodValue;
}







// Likelihood::getCovariates()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

ArrayXd Likelihood::getCovariates()
{
    return covariates;
}






// Likelihood::getObservations()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

ArrayXd Likelihood::getObservations()
{
    return observations;
}







// Likelihood::getUncertainties()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

ArrayXd Likelihood::getUncertainties()
{
    return uncertainties;
}






// Likelihood::getPredictions()
//
// PURPOSE:
//
// INPUT:
//
// OUTPUT:

ArrayXd Likelihood::getPredictions()
{
    return predictions;
}




