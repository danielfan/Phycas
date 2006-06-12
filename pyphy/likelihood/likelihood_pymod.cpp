#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle	// only needed on lewis.eeb.uconn.edu

#include "phycas/force_include.h"
#include "CipresCommlib/CipresDataMatrixHelper.h"
#include "phycas/rand/probability_distribution.hpp"
#include "pyphy/phylogeny/basic_tree.hpp"
#include "pyphy/likelihood/likelihood_models.hpp"
#include "pyphy/likelihood/tree_likelihood.hpp"
#include "pyphy/likelihood/tip_data.hpp"
#include "pyphy/likelihood/internal_data.hpp"
#include "pyphy/likelihood/mcmc_param.hpp"
#include "pyphy/likelihood/mcmc_chain_manager.hpp"
#include "pyphy/likelihood/topo_prior_calculator.hpp"
#include "pyphy/likelihood/larget_simon_move.hpp"
#include "pyphy/likelihood/ncat_move.hpp"
#include "pyphy/likelihood/bush_move.hpp"
#include "pyphy/likelihood/edge_move.hpp"
#include "pyphy/likelihood/sim_data.hpp"
#include "pyphy/likelihood/q_matrix.hpp"
#include "pyphy/likelihood/xlikelihood.hpp"

void more_likelihood_pymod();

#include <boost/python.hpp>
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

		int startTreeViewer(TreeShPtr t, std::string msg) const 
			{ 
			return call_method<int,TreeShPtr,std::string>(m_self, "startTreeViewer", t, msg); 
			}

	private:
		PyObject * const m_self;
	};

BOOST_PYTHON_MODULE(_LikelihoodBase)
{
	// these lines required by num_util
	import_array();
	numeric::array::set_module_and_type("numarray", "NDArray");

	class_<AdHocDensity, boost::noncopyable, boost::shared_ptr<AdHocDensity> >("AdHocDensityBase", no_init)
		;
	class_<phycas::MCMCChainManager, boost::noncopyable, boost::shared_ptr<phycas::MCMCChainManager> >("MCMCChainManagerBase")
		.def("finalize", &MCMCChainManager::finalize)
		.def("recalcEdgeLenPriors", &MCMCChainManager::recalcEdgeLenPriors)
		.def("calcJointLnPrior", &MCMCChainManager::calcJointLnPrior)
		.def("addMove", &MCMCChainManager::addMove)
		.def("addModelParam", &MCMCChainManager::addModelParam)
		.def("addEdgeLenParam", &MCMCChainManager::addEdgeLenParam)
		.def("addEdgeLenHyperparam", &MCMCChainManager::addEdgeLenHyperparam)
		.def("getLastLnLike", &MCMCChainManager::getLastLnLike)
		.def("getLastLnPrior", &MCMCChainManager::getLastLnPrior)
        .def("getAllUpdaters", &MCMCChainManager::getAllUpdaters, return_value_policy<copy_const_reference>())
        .def("getMoves", &MCMCChainManager::getMoves, return_value_policy<copy_const_reference>())
        .def("getModelParams", &MCMCChainManager::getModelParams, return_value_policy<copy_const_reference>())
        .def("getEdgeLenParams", &MCMCChainManager::getEdgeLenParams, return_value_policy<copy_const_reference>())
        .def("getEdgeLenHyperparam", &MCMCChainManager::getEdgeLenHyperparam)
		.def("addMCMCUpdaters", &MCMCChainManager::addMCMCUpdaters)
		.def("clear", &MCMCChainManager::clear)
		.def("refreshLastLnLike", &MCMCChainManager::refreshLastLnLike)
		.def("refreshLastLnPrior", &MCMCChainManager::refreshLastLnPrior)
		;
	class_<std::vector<MCMCUpdaterShPtr> >("paramVec", no_init)
		.def("__iter__",  iterator<std::vector<MCMCUpdaterShPtr> >())
		;
   	class_<phycas::MCMCUpdater, bases<AdHocDensity>, boost::noncopyable, boost::shared_ptr<phycas::MCMCUpdater> >("MCMCUpdaterBase", no_init)
		.def("getLnPrior", &MCMCUpdater::getLnPrior)
		.def("getLnLike", &MCMCUpdater::getLnLike)
		.def("getCurrValue", &MCMCUpdater::getCurrValue)
		.def("getName", &MCMCUpdater::getName, return_value_policy<copy_const_reference>())
		.def("getPriorDescr", &MCMCUpdater::getPriorDescr)
		.def("getWeight", &MCMCUpdater::getWeight)
		.def("setName", &MCMCUpdater::setName)
		.def("setWeight", &MCMCUpdater::setWeight)
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
		;
	class_<TreeLikelihood, TreeLikelihoodWrapper, boost::noncopyable>("TreeLikelihoodBase", init<boost::shared_ptr<Model> >())
		.def("copyDataFromDiscreteMatrix", &TreeLikelihood::copyDataFromDiscreteMatrix)
		.def("copyDataFromSimData", &TreeLikelihood::copyDataFromSimData)
		.def("prepareForSimulation", &TreeLikelihood::prepareForSimulation)
		.def("prepareForLikelihood", &TreeLikelihood::prepareForLikelihood)
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
		.def("setNoData", &TreeLikelihood::setNoData)
		.def("setHaveData", &TreeLikelihood::setHaveData)
		.def("getNPatterns", &TreeLikelihood::getNPatterns)
		.def("getLikelihoodRootNodeNum", &TreeLikelihood::getLikelihoodRootNodeNum)
		;
	class_<TipData, boost::noncopyable>("TipData", no_init)
		.def("parentalCLAValid", &TipData::parentalCLAValid)
		.def("parentalCLACached", &TipData::parentalCLACached)
		;
	class_<InternalData, boost::noncopyable>("InternalData", no_init)
		.def("filialCLAValid", &InternalData::filialCLAValid)
		.def("filialCLACached", &InternalData::filialCLACached)
		.def("parentalCLAValid", &InternalData::parentalCLAValid)
		.def("parentalCLACached", &InternalData::parentalCLACached)
		;
	class_<QMatrix, boost::noncopyable>("QMatrixBase")
		.def("getDimension", &QMatrix::getDimension)
		.def("setRelativeRates", &QMatrix::setRelativeRates)
		.def("setStateFreqs", &QMatrix::setStateFreqs)
		.def("getPMatrix", &QMatrix::getPMatrix)
		.def("getQMatrix", &QMatrix::getQMatrix)
		.def("getEigenVectors", &QMatrix::getEigenVectors)
		.def("getEigenValues", &QMatrix::getEigenValues)
		;

	// This function call necessary to avoid fatal error C1204: Compiler limit: internal structure overflow
	// in VC 7.1 (see http://www.boost.org/libs/python/doc/v2/faq.html#c1204)
	more_likelihood_pymod();

	register_exception_translator<XLikelihood>(&translateXLikelihood);
}
