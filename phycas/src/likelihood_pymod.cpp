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

#if defined(USING_NUMARRAY)
#	define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#endif

#include <boost/python.hpp>

//#include "phycas/force_include.h"
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/tip_data.hpp"
#include "phycas/src/internal_data.hpp"
#include "phycas/src/mcmc_param.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
//#include "phycas/src/topo_prior_calculator.hpp"
//#include "phycas/src/larget_simon_move.hpp"
//#include "phycas/src/ncat_move.hpp"
//#include "phycas/src/bush_move.hpp"
//#include "phycas/src/edge_move.hpp"
#include "phycas/src/sim_data.hpp"
#include "phycas/src/q_matrix.hpp"
#include "phycas/src/xlikelihood.hpp"

void model_pymod();
void updater_pymod();

using namespace boost::python;
using namespace phycas;

void translateXLikelihood(const XLikelihood &e)
	{
    // Use the Python 'C' API to set up an exception object
    PyErr_SetString(PyExc_Exception, e.what());
    }

// TreeLikelihoodWrapper is necessary in order to allow this (in C++ code)
//
// startTreeViewer(t);
//
// to be translated into a call of TreeLikelihoodDerived.startTreeViewer, where 
// TreeLikelihoodDerived is a Python class derived from TreeLikelihoodBase. 
// (TreeLikelihoodBase is what TreeLikelihood is called within Python, search for 
// TreeLikelihoodBase below.) For more background, go to 
// http://www.boost.org/libs/python/doc/index.html, click on the "Reference Manual"
// link, then see the example at the bottom of the page that results from clicking
// on the "call_method.hpp" link (under the category "Function Invocation and 
// Creation").
class TreeLikelihoodWrapper : public TreeLikelihood
	{
	public:
		TreeLikelihoodWrapper(PyObject * self, ModelShPtr m) : TreeLikelihood(m), m_self(self) 
			{
			}
		virtual ~TreeLikelihoodWrapper() {} //@POL may need to delete uMat here (see TreeLikelihood::~TreeLikelihood)

		int startTreeViewer(TreeShPtr t, std::string msg, unsigned site) const 
			{ 
			return call_method<int,TreeShPtr,std::string>(m_self, "startTreeViewer", t, msg, site); 
			}

	private:
		PyObject * const m_self;
	};

BOOST_PYTHON_MODULE(_LikelihoodBase)
{
#if defined(USING_NUMARRAY)
#error USING_NUMARRAY defined
	// these lines required by num_util
	import_array();
	numeric::array::set_module_and_type("numarray", "NDArray");
#endif

	class_<AdHocDensity, boost::noncopyable, boost::shared_ptr<AdHocDensity> >("AdHocDensityBase", no_init)
		;
	class_<phycas::MCMCChainManager, boost::noncopyable, boost::shared_ptr<phycas::MCMCChainManager> >("MCMCChainManagerBase")
		.def("finalize", &MCMCChainManager::finalize)
		.def("calcJointLnPrior", &MCMCChainManager::calcJointLnPrior)
		.def("addMove", &MCMCChainManager::addMove)
		.def("addModelParam", &MCMCChainManager::addModelParam)
		.def("addEdgeLenParam", &MCMCChainManager::addEdgeLenParam)
		.def("addEdgeLenHyperparam", &MCMCChainManager::addEdgeLenHyperparam)
		.def("setEdgeLenHyperparam", &MCMCChainManager::setEdgeLenHyperparam)
		.def("getLastLnLike", &MCMCChainManager::getLastLnLike)
		.def("getLastLnPrior", &MCMCChainManager::getLastLnPrior)
        .def("getAllUpdaters", &MCMCChainManager::getAllUpdaters, return_value_policy<copy_const_reference>())
        .def("getMoves", &MCMCChainManager::getMoves, return_value_policy<copy_const_reference>())
        .def("getModelParams", &MCMCChainManager::getModelParams, return_value_policy<copy_const_reference>())
        .def("getEdgeLenParams", &MCMCChainManager::getEdgeLenParams, return_value_policy<copy_const_reference>())
        .def("getEdgeLenHyperparams", &MCMCChainManager::getEdgeLenHyperparams)
		.def("addMCMCUpdaters", &MCMCChainManager::addMCMCUpdaters)
		.def("clear", &MCMCChainManager::clear)
		.def("refreshLastLnLike", &MCMCChainManager::refreshLastLnLike)
		.def("refreshLastLnPrior", &MCMCChainManager::refreshLastLnPrior)
		;
	class_<std::vector<MCMCUpdaterShPtr> >("paramVec", no_init)
		.def("__iter__",  iterator<std::vector<MCMCUpdaterShPtr> >())
		;
	class_<phycas::SimData, boost::noncopyable, boost::shared_ptr<phycas::SimData> >("SimDataBase")
		.def("clear", &phycas::SimData::clear)
		.def("zeroCounts", &phycas::SimData::zeroCounts)
		.def("createMapleTuples", &phycas::SimData::createMapleTuples)
		.def("appendCountsToFile", &phycas::SimData::appendCountsToFile)
		.def("getPatterns", &phycas::SimData::getPatterns)
		.def("getNUniquePatterns", &phycas::SimData::getNUniquePatterns)
		.def("resetPatternLength", &phycas::SimData::resetPatternLength)
		.def("wipePattern", &phycas::SimData::wipePattern)
		.def("setState", &phycas::SimData::setState)
		.def("insertPattern", &phycas::SimData::insertPattern)
		.def("saveToNexusFile", &phycas::SimData::saveToNexusFile)
		.def("patternTable", &phycas::SimData::patternTable)
		.def("divideBy", &phycas::SimData::divideBy)
		.def("addDataTo", &phycas::SimData::addDataTo)
		.def("calct", &phycas::SimData::calct)
		.def("multBy", &phycas::SimData::multBy)
		.def("debugAppendCountsToFile", &phycas::SimData::debugAppendCountsToFile)
		.def("calctBinned", &phycas::SimData::calctBinned)
		.def("buildBinVector", &phycas::SimData::buildBinVector)
		.def("getTotalCount", &phycas::SimData::getTotalCount)
        .def("getBinnedCounts", &phycas::SimData::getBinnedCounts)
		;
	class_<TreeLikelihood, TreeLikelihoodWrapper, boost::noncopyable>("TreeLikelihoodBase", init<boost::shared_ptr<Model> >())
        .def("calcLogLikeAtSubstitutionSaturation", &TreeLikelihood::calcLogLikeAtSubstitutionSaturation)
		.def("getPatternCounts", &TreeLikelihood::getPatternCounts, return_value_policy<copy_const_reference>())
		.def("getCharIndexToPatternIndex", &TreeLikelihood::getCharIndexToPatternIndex, return_value_policy<copy_const_reference>())
		.def("getSiteLikelihoods", &TreeLikelihood::getSiteLikelihoods, return_value_policy<copy_const_reference>())
		.def("getSiteUF", &TreeLikelihood::getSiteUF, return_value_policy<copy_const_reference>())
		.def("storingSiteLikelihoods", &TreeLikelihood::storingSiteLikelihoods)
		.def("storeSiteLikelihoods", &TreeLikelihood::storeSiteLikelihoods)
		.def("copyDataFromDiscreteMatrix", &TreeLikelihood::copyDataFromDiscreteMatrix)
		.def("copyDataFromSimData", &TreeLikelihood::copyDataFromSimData)
		.def("prepareForSimulation", &TreeLikelihood::prepareForSimulation)
		.def("prepareForLikelihood", &TreeLikelihood::prepareForLikelihood)
		.def("addOrphanTip", &TreeLikelihood::addOrphanTip)
		.def("addDecoratedInternalNode", &TreeLikelihood::addDecoratedInternalNode)
		.def("replaceModel", &TreeLikelihood::replaceModel)
		.def("invalidateAwayFromNode", &TreeLikelihood::invalidateAwayFromNode)
		.def("calcLnLFromNode", &TreeLikelihood::calcLnLFromNode)
		.def("calcLnL", &TreeLikelihood::calcLnL)
		.def("simulateFirst", &TreeLikelihood::simulateFirst)
		.def("simulate", &TreeLikelihood::simulate)
		.def("listPatterns", &TreeLikelihood::listPatterns)
		.def("resetNEvals", &TreeLikelihood::resetNEvals)
		.def("getNEvals", &TreeLikelihood::getNEvals)
		.def("addDataTo", &TreeLikelihood::addDataTo)
		.def("recalcRelativeRates", &TreeLikelihood::recalcRelativeRates)
		.def("getCategoryLowerBoundaries", &TreeLikelihood::getCategoryLowerBoundaries)
		.def("getRateMeans", &TreeLikelihood::getRateMeans, return_value_policy<copy_const_reference>())
		.def("getRateProbs", &TreeLikelihood::getRateProbs, return_value_policy<copy_const_reference>())
		.def("getListOfAllMissingSites", &TreeLikelihood::getListOfAllMissingSites, return_value_policy<copy_const_reference>())
		.def("setNoData", &TreeLikelihood::setNoData)
		.def("setHaveData", &TreeLikelihood::setHaveData)
		.def("getNPatterns", &TreeLikelihood::getNPatterns)
		.def("getLikelihoodRootNodeNum", &TreeLikelihood::getLikelihoodRootNodeNum)
		.def("useUnimap", &TreeLikelihood::useUnimap)
		.def("isUsingUnimap", &TreeLikelihood::isUsingUnimap)
		.def("fullRemapping", &TreeLikelihood::fullRemapping)
		.def("setUFNumEdges", &TreeLikelihood::setUFNumEdges)
		.def("bytesPerCLA", &TreeLikelihood::bytesPerCLA)
		.def("numCLAsCreated", &TreeLikelihood::numCLAsCreated)
		.def("numCLAsStored", &TreeLikelihood::numCLAsStored)
		.def("setLot", &TreeLikelihood::setLot)
		//.def("setDebug", &TreeLikelihood::setDebug)
		.def("debugCheckForUncachedCLAs", &TreeLikelihood::debugCheckForUncachedCLAs)
		;
	class_<TipData, boost::noncopyable>("TipData", no_init)
		.def("parentalCLAValid", &TipData::parentalCLAValid)
		.def("parentalCLACached", &TipData::parentalCLACached)
		.def("getNumUnivents", &TipData::getNumUnivents)
		.def("getUniventStates", &TipData::getUniventStates)
		.def("getUniventTimes", &TipData::getUniventTimes)
		;
	class_<InternalData, boost::noncopyable>("InternalData", no_init)
		.def("filialCLAValid", &InternalData::filialCLAValid)
		.def("filialCLACached", &InternalData::filialCLACached)
		.def("parentalCLAValid", &InternalData::parentalCLAValid)
		.def("parentalCLACached", &InternalData::parentalCLACached)
		.def("getNumUnivents", &InternalData::getNumUnivents)
		.def("getUniventStates", &InternalData::getUniventStates)
		.def("getUniventTimes", &InternalData::getUniventTimes)
		;
#if 0
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

	// This function call necessary to avoid fatal error C1204: Compiler limit: internal structure overflow
	// in VC 7.1 (see http://www.boost.org/libs/python/doc/v2/faq.html#c1204)
	model_pymod();
	updater_pymod();

	register_exception_translator<XLikelihood>(&translateXLikelihood);
}
