/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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
#include "phypy/src/cipres/CipresDataMatrixHelper.h"
#include "phypy/src/tree_likelihood.hpp"
#include "phypy/src/larget_simon_move.hpp"
#include "phypy/src/ncat_move.hpp"
#include "phypy/src/bush_move.hpp"
#include "phypy/src/edge_move.hpp"
#include "phypy/src/topo_prior_calculator.hpp"

#include <boost/python.hpp>
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
		;
	class_<phycas::KappaParam, bases<phycas::MCMCUpdater, AdHocDensity>, 
		boost::noncopyable, boost::shared_ptr<phycas::KappaParam> >("KappaParam")
		;
	class_<phycas::GTRRateParam, bases<phycas::MCMCUpdater, AdHocDensity>, 
		boost::noncopyable, boost::shared_ptr<phycas::GTRRateParam> >("GTRRateParam", init<unsigned>())
		;
	class_<phycas::BaseFreqParam, bases<phycas::MCMCUpdater, AdHocDensity>, 
		boost::noncopyable, boost::shared_ptr<phycas::BaseFreqParam> >("BaseFreqParam", init<unsigned>()) 
		;
	class_<phycas::HyperPriorParam, bases<phycas::MCMCUpdater, AdHocDensity>, 
		boost::noncopyable, boost::shared_ptr<phycas::HyperPriorParam> >("HyperPriorParam") 
		;
	class_<phycas::EdgeLenParam, bases<phycas::MCMCUpdater, AdHocDensity>, 
		boost::noncopyable, boost::shared_ptr<phycas::EdgeLenParam> >("EdgeLenParam", init<TreeNode *>()) 
		;
	class_<phycas::LargetSimonMove, bases<phycas::MCMCUpdater>, 
		boost::noncopyable, boost::shared_ptr<phycas::LargetSimonMove> >("LargetSimonMove") 
		.def("update", &phycas::LargetSimonMove::update)
		.def("setLambda", &phycas::LargetSimonMove::setLambda)
		.def("getLambda", &phycas::LargetSimonMove::getLambda)
		.def("topologyChanged", &phycas::LargetSimonMove::topologyChanged)
		.def("viewProposedMove", &phycas::LargetSimonMove::viewProposedMove)
		;
	class_<TopoPriorCalculator, boost::noncopyable, 
		boost::shared_ptr<phycas::TopoPriorCalculator> >("TopoPriorCalculatorBase")
		.def("setNTax", &TopoPriorCalculator::SetNTax)
		.def("getNTax", &TopoPriorCalculator::GetNTax)
		.def("chooseRooted", &TopoPriorCalculator::ChooseRooted)
		.def("chooseUnrooted", &TopoPriorCalculator::ChooseUnrooted)
		.def("getCount", &TopoPriorCalculator::GetCount)
		.def("getSaturatedCount", &TopoPriorCalculator::GetSaturatedCount)
		.def("getTotalCount", &TopoPriorCalculator::GetTotalCount)
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
		.def("viewProposedMove", &phycas::BushMove::viewProposedMove)
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
