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

#include <boost/python.hpp>

//#include "phycas/force_include.h"
#include "phycas/src/cipres/CipresDataMatrixHelper.h"
#include "phycas/src/basic_tree_node.hpp"
#include "phycas/src/mapping_move.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/larget_simon_move.hpp"
#include "phycas/src/unimap_nni_move.hpp"
#include "phycas/src/tree_scaler_move.hpp"
#include "phycas/src/ncat_move.hpp"
#include "phycas/src/bush_move.hpp"
#include "phycas/src/samc_move.hpp"
#include "phycas/src/edge_move.hpp"
#include "phycas/src/topo_prior_calculator.hpp"

using namespace boost::python;
using namespace phycas;

void updater_pymod()
	{
   	class_<phycas::MCMCUpdater, bases<AdHocDensity>, boost::noncopyable, boost::shared_ptr<phycas::MCMCUpdater> >("MCMCUpdaterBase", no_init)
		.def("getLnPrior", &MCMCUpdater::getLnPrior)
		.def("getLnLike", &MCMCUpdater::getLnLike)
		.def("getCurrValue", &MCMCUpdater::getCurrValue)
		.def("getName", &MCMCUpdater::getName, return_value_policy<copy_const_reference>())
		.def("getPriorDescr", &MCMCUpdater::getPriorDescr)
		.def("getWeight", &MCMCUpdater::getWeight)
		.def("setName", &MCMCUpdater::setName)
		.def("setWeight", &MCMCUpdater::setWeight)
		.def("setStartingValue", &MCMCUpdater::setStartingValue)
		.def("setTree", &MCMCUpdater::setTree)
		.def("setLot", &MCMCUpdater::setLot)
		.def("setTreeLikelihood", &MCMCUpdater::setTreeLikelihood)
		.def("setModel", &MCMCUpdater::setModel)
		.def("setChainManager", &MCMCUpdater::setChainManager)
		.def("hasSliceSampler", &MCMCUpdater::hasSliceSampler)
		.def("isParameter", &MCMCUpdater::isParameter)
		.def("isMasterParameter", &MCMCUpdater::isMasterParameter)
		.def("isHyperParameter", &MCMCUpdater::isHyperParameter)
		.def("isMove", &MCMCUpdater::isMove)
		.def("getSliceSampler", &MCMCUpdater::getSliceSampler)
		.def("update", &MCMCUpdater::update)
		.def("isFixed", &MCMCUpdater::isFixed)
		.def("fixParameter", &MCMCUpdater::fixParameter)
		.def("freeParameter", &MCMCUpdater::freeParameter)
		.def("getDebugInfo", &MCMCUpdater::getDebugInfo)
		.def("setSaveDebugInfo", &MCMCUpdater::setSaveDebugInfo)
		.def("setPower", &MCMCUpdater::setPower)
		.def("getPower", &MCMCUpdater::getPower)
		.def("setStandardHeating", &MCMCUpdater::setStandardHeating)
		.def("setLikelihoodHeating", &MCMCUpdater::setLikelihoodHeating)
		.def("isStandardHeating", &MCMCUpdater::isStandardHeating)
		.def("isLikelihoodHeating", &MCMCUpdater::isLikelihoodHeating)
		.def("isNoHeating", &MCMCUpdater::isNoHeating)
		;
	class_<phycas::KappaParam, bases<phycas::MCMCUpdater, AdHocDensity>, 
		boost::noncopyable, boost::shared_ptr<phycas::KappaParam> >("KappaParam")
		;
	class_<phycas::GTRRateParam, bases<phycas::MCMCUpdater, AdHocDensity>, 
		boost::noncopyable, boost::shared_ptr<phycas::GTRRateParam> >("GTRRateParam", init<unsigned>())
		;
	class_<phycas::StateFreqParam, bases<phycas::MCMCUpdater, AdHocDensity>, 
		boost::noncopyable, boost::shared_ptr<phycas::StateFreqParam> >("StateFreqParam", init<unsigned>()) 
		;
	class_<phycas::HyperPriorParam, bases<phycas::MCMCUpdater, AdHocDensity>, 
		boost::noncopyable, boost::shared_ptr<phycas::HyperPriorParam> >("HyperPriorParam") 
		;
	class_<phycas::LargetSimonMove, bases<phycas::MCMCUpdater>, 
		boost::noncopyable, boost::shared_ptr<phycas::LargetSimonMove> >("LargetSimonMove") 
		.def("update", &phycas::LargetSimonMove::update)
		.def("setLambda", &phycas::LargetSimonMove::setLambda)
		.def("getLambda", &phycas::LargetSimonMove::getLambda)
		.def("topologyChanged", &phycas::LargetSimonMove::topologyChanged)
		;
	class_<phycas::UnimapNNIMove, bases<phycas::MCMCUpdater>, 
		boost::noncopyable, boost::shared_ptr<phycas::UnimapNNIMove> >("UnimapNNIMove") 
		.def("update", &phycas::UnimapNNIMove::update)
		;
	class_<phycas::TreeScalerMove, bases<phycas::MCMCUpdater>, 
		boost::noncopyable, boost::shared_ptr<phycas::TreeScalerMove> >("TreeScalerMove") 
		.def("update", &phycas::TreeScalerMove::update)
		;
	class_<phycas::MappingMove, bases<phycas::MCMCUpdater>, 
		boost::noncopyable, boost::shared_ptr<phycas::MappingMove> >("MappingMove") 
		.def("update", &phycas::MappingMove::update)
		;
	class_<TopoPriorCalculator, boost::noncopyable, 
		boost::shared_ptr<phycas::TopoPriorCalculator> >("TopoPriorCalculatorBase")
		.def("setNTax", &TopoPriorCalculator::SetNTax)
		.def("getNTax", &TopoPriorCalculator::GetNTax)
		.def("chooseRooted", &TopoPriorCalculator::ChooseRooted)
		.def("chooseUnrooted", &TopoPriorCalculator::ChooseUnrooted)
        .def("getLnCount", &TopoPriorCalculator::GetLnCount)
		.def("getLnSaturatedCount", &TopoPriorCalculator::GetLnSaturatedCount)
		.def("getLnTotalCount", &TopoPriorCalculator::GetLnTotalCount)
		.def("getNFactorsVect", &TopoPriorCalculator::GetNFactorsVect)
		.def("getLnCounts", &TopoPriorCalculator::GetLnCounts)
		.def("setLnScalingFactor", &TopoPriorCalculator::SetLnScalingFactor)
		.def("getLnScalingFactor", &TopoPriorCalculator::GetLnScalingFactor)
		.def("getCountsVect", &TopoPriorCalculator::GetCountsVect)
		.def("chooseResolutionClassPrior", &TopoPriorCalculator::ChooseResolutionClassPrior)
		.def("choosePolytomyPrior", &TopoPriorCalculator::ChoosePolytomyPrior)
		.def("setC", &TopoPriorCalculator::SetC)
		.def("getC", &TopoPriorCalculator::GetC)
		.def("getLnTopologyPrior", &TopoPriorCalculator::GetLnTopologyPrior)
		.def("getLnNormalizedTopologyPrior", &TopoPriorCalculator::GetLnNormalizedTopologyPrior)
		.def("GetLnNormConstant", &TopoPriorCalculator::GetLnNormConstant)
		.def("getTopoPriorVect", &TopoPriorCalculator::GetTopoPriorVect)
		.def("isResolutionClassPrior", &TopoPriorCalculator::IsResolutionClassPrior)
		.def("isPolytomyPrior", &TopoPriorCalculator::IsPolytomyPrior)
		.def("isRooted", &TopoPriorCalculator::IsRooted)
		.def("isUnrooted", &TopoPriorCalculator::IsUnrooted)
		.def("getRealizedResClassPriorsVect", &TopoPriorCalculator::GetRealizedResClassPriorsVect)
		;
	class_<phycas::BushMove, bases<phycas::MCMCUpdater>, 
		boost::noncopyable, boost::shared_ptr<phycas::BushMove> >("BushMove") 
		.def("update", &phycas::BushMove::update)
		.def("addEdgeMoveProposed", &phycas::BushMove::addEdgeMoveProposed)
		.def("setEdgeLenDistMean", &phycas::BushMove::setEdgeLenDistMean)
		.def("finalize", &phycas::BushMove::finalize)
		.def("getTopoPriorCalculator", &phycas::BushMove::getTopoPriorCalculator)
		;
	class_<phycas::SamcMove, bases<phycas::MCMCUpdater>, 
		boost::noncopyable, boost::shared_ptr<phycas::SamcMove> >("SamcMove", init<ProbDistShPtr>()) 
		.def("extrapolate", &phycas::SamcMove::extrapolate)
		.def("project", &phycas::SamcMove::project)
		.def("update", &phycas::SamcMove::update)
		//.def("finalize", &phycas::SamcMove::finalize)
		.def("setDistanceMatrixRow", &phycas::SamcMove::setDistanceMatrixRow)
		.def("setTemperature", &phycas::SamcMove::setTemperature)
		.def("getTemperature", &phycas::SamcMove::getTemperature)
		.def("setNTax", &phycas::SamcMove::setNTax)
		.def("getNTax", &phycas::SamcMove::getNTax)
		;
	class_<phycas::NCatMove, bases<phycas::MCMCUpdater>, 
		boost::noncopyable, boost::shared_ptr<phycas::NCatMove> >("NCatMove") 
		.def("update", &phycas::NCatMove::update)
		.def("addCatMoveProposed", &phycas::NCatMove::addCatMoveProposed)
		.def("getPhi", &phycas::NCatMove::getPhi)
		.def("setPhi", &phycas::NCatMove::setPhi)
		.def("getLambda", &phycas::NCatMove::getLambda)
		.def("setLambda", &phycas::NCatMove::setLambda)
		.def("getNCatMax", &phycas::NCatMove::getNCatMax)
		.def("setNCatMax", &phycas::NCatMove::setNCatMax)
		.def("getS", &phycas::NCatMove::getS)
		.def("setS", &phycas::NCatMove::setS)
		.def("getL", &phycas::NCatMove::getL)
		.def("setL", &phycas::NCatMove::setL)
		.def("getCatProbPrior", &phycas::NCatMove::getCatProbPrior)
		.def("setCatProbPrior", &phycas::NCatMove::setCatProbPrior)
		;
	class_<phycas::EdgeMove, bases<phycas::MCMCUpdater>, 
		boost::noncopyable, boost::shared_ptr<phycas::EdgeMove> >("EdgeMove") 
		.def("update", &phycas::EdgeMove::update)
		.def("setLambda", &phycas::EdgeMove::setLambda)
		.def("getLambda", &phycas::EdgeMove::getLambda)
		;
	}
