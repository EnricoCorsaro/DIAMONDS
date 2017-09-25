//
// Created by mamu on 9/20/17.
//
//pybind includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/chrono.h>
//Diamonds includes
#include "UniformPrior.h"
#include "Model.h"
#include "Clusterer.h"
#include "Metric.h"
#include "ExponentialLikelihood.h"
#include "Ellipsoid.h"
#include "EuclideanMetric.h"
#include "FerozReducer.h"
#include "GridUniformPrior.h"
#include "KmeansClusterer.h"
#include "Likelihood.h"
#include "LivePointsReducer.h"
#include "NestedSampler.h"
#include "MeanNormalLikelihood.h"
#include "Metric.h"
#include "MultiEllipsoidSampler.h"
#include "NormalLikelihood.h"
#include "NormalPrior.h"
#include "PowerlawReducer.h"
#include "Prior.h"
#include "Results.h"
#include "SuperGaussianPrior.h"
#include "ZeroClusterer.h"
#include "ZeroModel.h"
#include "ZeroPrior.h"
#include "ZeroSampler.h"
//Python_wrapper includes
#include "include/PyModel.h"
#include "include/PyClusterer.h"
#include "include/PyEuclideanMetric.h"
#include "include/PyFerozReducer.h"
#include "include/PyGridUniformPrior.h"
#include "include/PyLikelihood.h"
#include "include/PyLivePointsReducer.h"
#include "include/PyMeanNormalLikelihood.h"
#include "include/PyMetric.h"
#include "include/PyMultiEllispoidSampler.h"
#include "include/PyNestedSampler.h"
#include "include/PyNormalLikelihood.h"
#include "include/PyNormalPrior.h"
#include "include/PyPowerlawReducer.h"
#include "include/PyPrior.h"
#include "include/PySuperGaussianPrior.h"
#include "include/PyUniformPrior.h"
#include "include/PyZeroModel.h"
#include "include/PyZeroPrior.h"
#include "include/PyZeroSampler.h"

namespace py = pybind11;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;
typedef Eigen::Ref<Eigen::ArrayXXd> RefArrayXXd;

PYBIND11_MODULE(pyDiamonds,m)
{
    py::class_<Clusterer,PyClusterer<>> (m,"Clusterer")
        .def(py::init_alias<Metric&>())
        .def("cluster",&Clusterer::cluster);

    py::class_<Ellipsoid>(m,"Ellipsoid")
        .def(py::init<RefArrayXXd,const double>())
        .def("resetEnlargementFraction",&Ellipsoid::resetEnlargementFraction)
        .def("overlapsWith",&Ellipsoid::overlapsWith)
        .def("containsPoint",&Ellipsoid::containsPoint)
        .def("drawPoint",&Ellipsoid::drawPoint)
        .def("getCenterCoordinates",&Ellipsoid::getCenterCoordinates)
        .def("getEigenvalues",&Ellipsoid::getEigenvalues)
        .def("getSample",&Ellipsoid::getSample)
        .def("getCovarianceMatrix",&Ellipsoid::getCovarianceMatrix)
        .def("getEigenvectors",&Ellipsoid::getEigenvectors)
        .def("getSampleSize",&Ellipsoid::getSampleSize)
        .def("getHyperVolume",&Ellipsoid::getHyperVolume)
        .def("getEnlargementFraction",&Ellipsoid::getEnlargementFraction);

    py::class_<Metric,PyMetric<>>(m,"Metric")
            .def(py::init_alias<>())
            .def("distance",&Metric::distance);

    py::class_<Likelihood,PyLikelihood>(m,"Likelihood")
            .def(py::init_alias<const RefArrayXd,Model&>())
            .def("logValue",&Likelihood::logValue);

    py::class_<Model,PyModel> (m,"Model")
            .def(py::init_alias<const RefArrayXd>())
            .def("getCovariates",&Model::getCovariates)
            .def("getNparameters",&Model::getNparameters)
            .def("predict",&Model::predict);

    py::class_<Prior,PyPrior>(m,"Prior")
            .def(py::init<const int>())
            .def("getNdimensions",&Prior::getNdimensions)
            .def("logDensity",&Prior::logDensity)
            .def("drawnPointIsAccepted",&Prior::drawnPointIsAccepted)
            .def("draw",&Prior::draw)
            .def("drawWithConstraint",&Prior::drawWithConstraint)
            .def("writeHyperParametersToFile",&Prior::writeHyperParametersToFile);

    py::class_<LivePointsReducer,PyLivePointsReducer>(m,"LivePointsReducer")
            .def(py::init_alias<NestedSampler&>())
            .def("findIndicesOfLivePointsToRemove",&LivePointsReducer::findIndicesOfLivePointsToRemove)
            .def("getNlivePointsToRemove",&LivePointsReducer::getNlivePointsToRemove)
            .def("updateNlivePoints",&LivePointsReducer::updateNlivePoints);

    py::class_<NestedSampler,PyNestedSampler>(m,"NestedSampler")
            .def(py::init<const bool, const int, const int, vector<Prior*>,Likelihood&,Metric&,Clusterer&>())
            .def("run",&NestedSampler::run)
            .def("drawWithConstraint",&NestedSampler::drawWithConstraint)
            .def("getNiterations",&NestedSampler::getNiterations)
            .def("getNdimensions",&NestedSampler::getNdimensions)
            .def("getNlivePoints",&NestedSampler::getNlivePoints)
            .def("getInitialNlivePoints",&NestedSampler::getInitialNlivePoints)
            .def("getMinNlivePoints",&NestedSampler::getMinNlivePoints)
            .def("getLogCumulatedPriorMass",&NestedSampler::getLogCumulatedPriorMass)
            .def("getLogRemainingPriorMass",&NestedSampler::getLogRemainingPriorMass)
            .def("getRatioOfRemainderToCurrentEvidence",&NestedSampler::getRatioOfRemainderToCurrentEvidence)
            .def("getLogMaxLikelihoodOfLivePoints",&NestedSampler::getLogMaxLikelihoodOfLivePoints)
            .def("getComputationalTime",&NestedSampler::getComputationalTime)
            .def("getTerminationFactor",&NestedSampler::getTerminationFactor)
            .def("getNlivePointsPerIteration",&NestedSampler::getNlivePointsPerIteration)
            .def("getNestedSample",&NestedSampler::getNestedSample)
            .def("getLogLikelihood",&NestedSampler::getLogLikelihood)
            .def("setLogEvidence",&NestedSampler::setLogEvidence)
            .def("getLogEvidence",&NestedSampler::getLogEvidence)
            .def("setLogEvidenceError",&NestedSampler::setLogEvidenceError)
            .def("getLogEvidenceError",&NestedSampler::getLogEvidenceError)
            .def("setInformationGain",&NestedSampler::setInformationGain)
            .def("getInformationGain",&NestedSampler::getInformationGain)
            .def("setPosteriorSample",&NestedSampler::setPosteriorSample)
            .def("getPosteriorSample",&NestedSampler::getPosteriorSample)
            .def("setLogLikelihoodOfPosteriorSample",&NestedSampler::setLogLikelihoodOfPosteriorSample)
            .def("getLogLikelihoodOfPosteriorSample",&NestedSampler::getLogLikelihoodOfPosteriorSample)
            .def("setLogWeightOfPosteriorSample",&NestedSampler::setLogWeightOfPosteriorSample)
            .def("getLogWeightOfPosteriorSample",&NestedSampler::getLogWeightOfPosteriorSample)
            .def("setOutputPathPrefix",&NestedSampler::setOutputPathPrefix)
            .def("getOutputPathPrefix",&NestedSampler::getOutputPathPrefix)
            .def("verifySamplerStatus",&NestedSamplerPublicist::verifySamplerStatus);

    py::class_<EuclideanMetric,PyMetric<EuclideanMetric>,Metric>(m,"EuclideanMetric")
        .def(py::init_alias<>())
        .def("distance",&EuclideanMetric::distance);

    py::class_<ExponentialLikelihood,Likelihood>(m,"ExponentialLikelihood")
        .def(py::init<const RefArrayXd,Model&>())
        .def("logValue",&ExponentialLikelihood::logValue);

    py::class_<FerozReducer,PyFerozReducer,LivePointsReducer>(m,"FerozReducer")
        .def(py::init_alias<NestedSampler&,const double>())
        .def("updateNlivePoints",&FerozReducer::updateNlivePoints);

    py::class_<GridUniformPrior,PyGridUniformPrior,Prior>(m,"GridUniformPrior")
        .def(py::init_alias<const RefArrayXd,const RefArrayXd,const RefArrayXd, const RefArrayXd>())
        .def("getStartingCoordinate",&GridUniformPrior::getStartingCoordinate)
        .def("getNgridPoints",&GridUniformPrior::getNgridPoints)
        .def("getSeparation",&GridUniformPrior::getSeparation)
        .def("getTolerance",&GridUniformPrior::getTolerance)
        .def("logDensity",&GridUniformPrior::logDensity)
        .def("drawnPointIsAccepted",&GridUniformPrior::drawnPointIsAccepted)
        .def("draw",&GridUniformPrior::draw)
        .def("drawWithConstraint",&GridUniformPrior::drawWithConstraint)
        .def("writeHyperParametersToFile",&GridUniformPrior::writeHyperParametersToFile);

    py::class_<KmeansClusterer,PyClusterer<KmeansClusterer>,Clusterer>(m,"KmeansClusterer")
        .def(py::init_alias<Metric&,unsigned int, unsigned int, unsigned int, double>())
        .def("cluster",&KmeansClusterer::cluster);

    /*
    py::class_<MeanNormalLikelihood,PyMeanNormalLikelihood,Likelihood>(m,"MeanNormalLikelihood")
        .def(py::init_alias<const RefArrayXd,const RefArrayXd,Model&>())
        .def("getUncertainties",&MeanNormalLikelihood::getUncertainties)
        .def("getNormalizedUncertainties",&MeanNormalLikelihood::getNormalizedUncertainties)
        .def("getWeights",&MeanNormalLikelihood::getWeights)
        .def("logValue",&MeanNormalLikelihood::logValue);
    */

    py::class_<MultiEllipsoidSampler,PyMultiEllipsoidSampler,NestedSampler>(m,"MultiEllipsoidSampler")
        .def(py::init<const bool,vector<Prior*>,Likelihood&,Metric&,Clusterer &,const int, const int, const double, const double>())
        .def("drawWithConstraint",&MultiEllipsoidSampler::drawWithConstraint)
        .def("verifySamplerStatus",&MultiEllipsoidSampler::verifySamplerStatus)
        .def("getEllipsoids",&MultiEllipsoidSampler::getEllipsoids)
        .def("getInitialEnlargementFraction",&MultiEllipsoidSampler::getInitialEnlargementFraction)
        .def("getShrinkingRate",&MultiEllipsoidSampler::getShrinkingRate);

    py::class_<NormalLikelihood,PyNormalLikelihood,Likelihood>(m,"NormalLikelihood")
        .def(py::init<const RefArrayXd,const RefArrayXd,Model&>())
        .def("getUncertainties",&NormalLikelihood::getUncertainties)
        .def("logValue",&NormalLikelihood::logValue);

    py::class_<NormalPrior,PyNormalPrior,Prior>(m,"NormalPrior")
        .def(py::init_alias<RefArrayXd const,RefArrayXd const>())
        .def("getMean",&NormalPrior::getMean)
        .def("getStandardDeviation",&NormalPrior::getStandardDeviation)
        .def("logDensity",&NormalPrior::logDensity)
        .def("drawnPointIsAccepted",&NormalPrior::drawnPointIsAccepted)
        .def("draw",&NormalPrior::draw)
        .def("drawWithConstraint",&NormalPrior::drawWithConstraint)
        .def("writeHyperParametersToFile",&NormalPrior::writeHyperParametersToFile);

    py::class_<PowerlawReducer,PyPowerlawReducer,LivePointsReducer>(m,"PowerlawReducer")
        .def(py::init<NestedSampler&,const double,const double,const double>())
        .def("updateNlivePoints",&PowerlawReducer::updateNlivePoints);


    py::class_<Results>(m,"Results")
        .def(py::init<NestedSampler&>())
        .def("writeParametersToFile",&Results::writeParametersToFile)
        .def("writeLogLikelihoodToFile",&Results::writeLogLikelihoodToFile)
        .def("writeLogWeightsToFile",&Results::writeLogWeightsToFile)
        .def("writeEvidenceInformationToFile",&Results::writeEvidenceInformationToFile)
        .def("writePosteriorProbabilityToFile",&Results::writePosteriorProbabilityToFile)
        .def("writeParametersSummaryToFile",&Results::writeParametersSummaryToFile)
        .def("writeObjectsIdentificationToFile",&Results::writeObjectsIdentificationToFile);

    py::class_<SuperGaussianPrior,PySuperGaussianPrior,Prior>(m,"SuperGaussianPrior")
        .def(py::init<const RefArrayXd, const RefArrayXd, const RefArrayXd>())
        .def("getCenter",&SuperGaussianPrior::getCenter)
        .def("getSigma",&SuperGaussianPrior::getSigma)
        .def("getWidthOfPlateau",&SuperGaussianPrior::getWidthOfPlateau)
        .def("logDensity",&SuperGaussianPrior::logDensity)
        .def("drawnPointIsAccepted",&SuperGaussianPrior::drawnPointIsAccepted)
        .def("draw",&SuperGaussianPrior::draw)
        .def("drawWithConstraint",&SuperGaussianPrior::drawWithConstraint)
        .def("writeHyperParametersToFile",&SuperGaussianPrior::writeHyperParametersToFile);

    py::class_<UniformPrior,PyUniformPrior,Prior>(m,"UniformPrior")
        .def(py::init<const RefArrayXd,const RefArrayXd >())
        .def("getMinima", &UniformPrior::getMinima)
        .def("getMaxima", &UniformPrior::getMaxima)
        .def("logDensity",&UniformPrior::logDensity)
        .def("drawnPointIsAccepted",&UniformPrior::drawnPointIsAccepted)
        .def("draw",&UniformPrior::draw)
        .def("drawWithConstraint",&UniformPrior::drawWithConstraint)
        .def("writeHyperParametersToFile",&UniformPrior::writeHyperParametersToFile);

    py::class_<ZeroClusterer,PyClusterer<ZeroClusterer>,Clusterer>(m,"ZeroClusterer")
        .def(py::init<Metric&>())
        .def("cluster",&ZeroClusterer::cluster);

    py::class_<ZeroModel,PyZeroModel,Model>(m,"ZeroModel")
        .def(py::init<const RefArrayXd>())
        .def("predict",&ZeroModel::predict);

    py::class_<ZeroPrior,PyZeroPrior,Prior>(m,"ZeroPrior")
        .def(py::init<const int>())
        .def("logDensity",&ZeroPrior::logDensity)
        .def("drawnPointIsAccepted",&ZeroPrior::drawnPointIsAccepted)
        .def("draw",&ZeroPrior::draw)
        .def("drawWithConstraint",&ZeroPrior::drawWithConstraint)
        .def("writeHyperParametersToFile",&ZeroPrior::writeHyperParametersToFile);
/*
    py::class_<ZeroSampler,PyZeroSampler>(m,"ZeroSampler")
        .def(py::init<const bool, const int, const int, vector<Prior*>,Likelihood&, Metric&, Clusterer&>())
        .def("drawWithConstraint",&ZeroSampler::drawWithConstraint)
        .def("verifySamplerStatus",&NestedSamplerPublicist::verifySamplerStatus);
*/

}
