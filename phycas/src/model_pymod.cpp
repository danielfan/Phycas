/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

//#include "phycas/force_include.h"
#include "phycas/src/likelihood_models.hpp"

#include <boost/python.hpp>
using namespace boost::python;
using namespace phycas;

void model_pymod()
	{
#if 1
	class_<QMatrix, boost::noncopyable>("QMatrixBase")
		.def("getDimension", &QMatrix::getDimension)
		.def("setRelativeRates", &QMatrix::setRelativeRates)
		.def("setStateFreqs", &QMatrix::setStateFreqs)
		.def("getPMatrix", &QMatrix::getPMatrix)
		.def("getQMatrix", &QMatrix::getQMatrix)
		.def("getEigenVectors", &QMatrix::getEigenVectors)
		.def("getEigenValues", &QMatrix::getEigenValues)
		;
#endif
	class_<phycas::Model, boost::noncopyable>("Model", no_init)
		.def("fixEdgeLenHyperprior", &phycas::Model::fixEdgeLenHyperprior)
		.def("freeEdgeLenHyperprior", &phycas::Model::freeEdgeLenHyperprior)
		.def("fixEdgeLengths", &phycas::Model::fixEdgeLengths)
		.def("freeEdgeLengths", &phycas::Model::freeEdgeLengths)
		.def("fixStateFreqs", &phycas::Model::fixStateFreqs)
		.def("freeStateFreqs", &phycas::Model::freeStateFreqs)
		.def("fixShape", &phycas::Model::fixShape)
		.def("freeShape", &phycas::Model::freeShape)
		.def("shapeFixed", &phycas::Model::shapeFixed)
		.def("pinvarFixed", &phycas::Model::pinvarFixed)
		.def("stateFreqsFixed", &phycas::Model::stateFreqsFixed)
		.def("edgeLenHyperParamFixed", &phycas::Model::edgeLenHyperParamFixed)
		.def("edgeLengthsFixed", &phycas::Model::edgeLengthsFixed)
		.def("setFlexRateUpperBound", &phycas::Model::setFlexRateUpperBound)
		.def("setNumFlexSpacers", &phycas::Model::setNumFlexSpacers)
		.def("setFlexModel", &phycas::Model::setFlexModel)
		.def("setNotFlexModel", &phycas::Model::setNotFlexModel)
		.def("fixFlexProbs", &phycas::Model::fixFlexProbs)
		.def("freeFlexProbs", &phycas::Model::freeFlexProbs)
		.def("fixFlexRates", &phycas::Model::fixFlexRates)
		.def("freeFlexRates", &phycas::Model::freeFlexRates)
		.def("setFLEXProbParamPrior", &phycas::Model::setFLEXProbParamPrior)
		.def("getFLEXProbParamPrior", &phycas::Model::getFLEXProbParamPrior)
		.def("setFlexRateUnnorm", &phycas::Model::setFlexRateUnnorm)
		.def("setFlexProbUnnorm", &phycas::Model::setFlexProbUnnorm)
		.def("setPinvarModel", &phycas::Model::setPinvarModel)
		.def("setNotPinvarModel", &phycas::Model::setNotPinvarModel)
		.def("fixPinvar", &phycas::Model::fixPinvar)
		.def("freePinvar", &phycas::Model::freePinvar)
		.def("getPinvarPrior", &phycas::Model::getPinvarPrior)
		.def("setPinvarPrior", &phycas::Model::setPinvarPrior)
		.def("getDiscreteGammaShapePrior", &phycas::Model::getDiscreteGammaShapePrior)
		.def("setDiscreteGammaShapePrior", &phycas::Model::setDiscreteGammaShapePrior)
		.def("getEdgeLenPrior", &phycas::Model::getEdgeLenPrior)
		.def("setEdgeLenPrior", &phycas::Model::setEdgeLenPrior)
		.def("getEdgeLenHyperPrior", &phycas::Model::getEdgeLenHyperPrior)
		.def("setEdgeLenHyperPrior", &phycas::Model::setEdgeLenHyperPrior)
		.def("getModelName", &phycas::Model::getModelName)
		.def("getPinvar", &phycas::Model::getPinvar)
		.def("setPinvar", &phycas::Model::setPinvar)
		.def("getShape", &phycas::Model::getShape)
		.def("setShape", &phycas::Model::setShape)
		.def("setPriorOnShapeInverse", &phycas::Model::setPriorOnShapeInverse)
		.def("getNGammaRates", &phycas::Model::getNGammaRates)
		.def("setNGammaRates", &phycas::Model::setNGammaRates)
		.def("getGammaRateProbs", &phycas::Model::getGammaRateProbs, return_value_policy<copy_const_reference>())
		.def("setAllGammaRateProbsEqual", &phycas::Model::setAllGammaRateProbsEqual)
		.def("getPMatrix", &phycas::Model::getPMatrix)
		;
	class_<phycas::JC, bases<phycas::Model> >("JCModelBase")
		.def("getModelName", &phycas::JC::getModelName)
		.def("getNStates", &phycas::JC::getNStates)
		.def("getStateFreqs", &phycas::JC::getStateFreqs, return_value_policy<copy_const_reference>())
		.def("setAllFreqsEqual", &phycas::JC::setAllFreqsEqual)
		.def("paramHeader", &phycas::JC::paramHeader)
		.def("paramReport", &phycas::JC::paramReport)
		;
	class_<phycas::HKY, bases<phycas::Model> >("HKYModelBase")
		.def("getModelName", &phycas::HKY::getModelName)
		.def("fixKappa", &phycas::HKY::fixKappa)
		.def("freeKappa", &phycas::HKY::freeKappa)
		.def("getKappa", &phycas::HKY::getKappa)
		.def("setKappa", &phycas::HKY::setKappa)
		.def("getKappaPrior", &phycas::HKY::getKappaPrior)
		.def("setKappaPrior", &phycas::HKY::setKappaPrior)
		.def("getBaseFreqParamPrior", &phycas::HKY::getBaseFreqParamPrior)
		.def("setBaseFreqParamPrior", &phycas::HKY::setBaseFreqParamPrior)
		.def("setKappaFromTRatio", &phycas::HKY::setKappaFromTRatio)
		.def("calcTRatio", &phycas::HKY::calcTRatio)
		.def("getNStates", &phycas::HKY::getNStates)
		.def("getStateFreqs", &phycas::HKY::getStateFreqs, return_value_policy<copy_const_reference>())
		.def("setAllFreqsEqual", &phycas::HKY::setAllFreqsEqual)
		.def("setNucleotideFreqs", &phycas::HKY::setNucleotideFreqs)
		.def("setStateFreqUnnorm", &phycas::HKY::setStateFreqUnnorm)
		.def("paramHeader", &phycas::HKY::paramHeader)
		.def("paramReport", &phycas::HKY::paramReport)
		;
	class_<phycas::GTR, bases<phycas::Model> >("GTRModelBase")
		.def("getModelName", &phycas::GTR::getModelName)
		//.def("createParameters", &phycas::GTR::createParameters)
		//.def("calcPMat", &phycas::GTR::calcPMat)
		.def("fixRelRates", &phycas::GTR::fixRelRates)
		.def("freeRelRates", &phycas::GTR::freeRelRates)
		.def("getRelRates", &phycas::GTR::getRelRates)
		.def("setRelRates", &phycas::GTR::setRelRates)
		.def("setRelRateUnnorm", &phycas::GTR::setRelRateUnnorm)
		.def("setRelRatePrior", &phycas::GTR::setRelRatePrior)
		.def("getRelRatePrior", &phycas::GTR::getRelRatePrior)
		.def("setNucleotideFreqs", &phycas::GTR::setNucleotideFreqs)
		.def("setAllFreqsEqual", &phycas::GTR::setAllFreqsEqual)
		.def("setStateFreqUnnorm", &phycas::GTR::setStateFreqUnnorm)
		.def("setBaseFreqParamPrior", &phycas::GTR::setBaseFreqParamPrior)
		.def("getBaseFreqParamPrior", &phycas::GTR::getBaseFreqParamPrior)
		.def("paramHeader", &phycas::GTR::paramHeader)
		.def("paramReport", &phycas::GTR::paramReport)
		.def("calcTRatio", &phycas::GTR::calcTRatio)
		;
	}
