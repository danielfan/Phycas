#ifndef PHO_MCMC_H
#define PHO_MCMC_H

//#define CHECK_WITH_PAUP

class NxsCommandManager;
class DistributionManager;
class LongOperationManager;
class PhoTaxaManager;
class PhoCharactersManager;
class ProbabilityDistribution;
class LargetSimonMove;
class BushMaster;
class ExponentialDistribution;
class KappaDistribution;
class EdgeMove;
class TreeScalerMove;
class Tree;
class SliceSampler;
class GibbsParameter;
class GibbsHyperParameter;

class Lot;
class SplitManager;
class TopologyManager;
class EdgeLenPriorCalculator;
class HKYAdHocEvaluator;
class MCMCMove;
#if defined (CORBA_CLIENT_PHYCAS)
	namespace cipresCORBA
		{
		class PutTreeWrapper;
		}
#endif

#include "phycas/misc/stopwatch.hpp"
#include "phycas/command/mcmc_settings.hpp"
#include "phycas/modules/mcmc/move_schedule.hpp"
#include "phycas/modules/mcmc/prior_summary.hpp"
#include "phycas/modules/mcmc/chain_state.hpp"
#include <boost/shared_ptr.hpp>
typedef boost::shared_ptr<Lot>						LotShPtr;
typedef boost::shared_ptr<SplitManager>				SplitManagerShPtr;
typedef boost::shared_ptr<TopologyManager>			TopologyManagerShPtr;
typedef boost::shared_ptr<EdgeLenPriorCalculator>	EdgeLenPriorCalculatorShPtr;
typedef boost::shared_ptr<HKYAdHocEvaluator>		HKYAdHocEvaluatorShPtr;
typedef boost::shared_ptr<GibbsParameter>			GibbsParameterShPtr;
typedef boost::shared_ptr<GibbsHyperParameter>		GibbsHyperParameterShPtr;
typedef boost::shared_ptr<ProbabilityDistribution>	ProbabilityDistributionShPtr;
#include "ncl/output/nxs_output_destination_description.hpp" // for storing the destination of the split summary

/*----------------------------------------------------------------------------------------------------------------------
|	Class that handles the MCMC command and orchestrates an MCMC analysis.
*/
class PhoMCMC
	{
	public:
		
		CmdResult 						HandleMCMC(MCMCSettings *s);
		void 							SetupMCMCCommand(NxsCommandManager *cmdMgr);

	private:

		void							PrepareSampler();
		void							ShowInfo();
		bool							StartingSampler();
		void							RunSampler();
		void							StoppingSampler();
		void							NextStep();
		void							GibbsSampling(bool ignore_data = false);
		void							PlotBestTree();
		void							GibbsSampleParam(GibbsParameter * param, bool ignoreData);
		void							SampleChainState(const ChainState  & ) const;
		const ChainState			  & RefreshInternalCopyOfChainState();
		void							ReportChainState(const ChainState & ) const;

		void							DebugCheckWithPAUP();

	private:
		StopWatch						stopwatch;				/**< keeps track of elapsed wall clock time needed for each run */
		LotShPtr						r;						/**< pseudorandom number generator */
		SplitManagerShPtr				splitMgr;				/**< keeps track of splits encountered in sampled trees */
		TopologyManagerShPtr			topoMgr;				/**< keeps track of topologies encountered in sampled trees */
		HKYAdHocEvaluatorShPtr			hky;					/**< pointer to object that handles likelihood calculations */

		MoveSchedule					moveSchedule;			///
		ChainState						lastRecordedChainState; ///
		PriorSummary					topoPrior;
		ProbabilityDistributionShPtr	edgelenPrior;			/**< prior distribution assumed for edge lengths */
		ProbabilityDistributionShPtr	freqParamPrior;			/**< prior distribution assumed for frequency parameters */
		ProbabilityDistributionShPtr	kappaPrior;				/**< prior distribution assumed for kappa */
		ProbabilityDistributionShPtr	rateVarPrior;			/**< prior distribution assumed for the variance in relative rates */
		ProbabilityDistributionShPtr	edgelenHyperPrior;		/**< hyperprior distribution assumed for the edge length distribution */
		ProbabilityDistributionShPtr	exponDistr;				/**< provides random deviates for edge length proposals */
		EdgeLenPriorCalculatorShPtr		edgeLenPriorCalc;		/**< computes prior for edge lengths given a probability distribution */

		GibbsParameterShPtr				freqAParam;				/**< updateable dirichlet parameter governing freq. A */
		GibbsParameterShPtr				freqCParam;				/**< updateable dirichlet parameter governing freq. C */
		GibbsParameterShPtr				freqGParam;				/**< updateable dirichlet parameter governing freq. G */
		GibbsParameterShPtr				freqTParam;				/**< updateable dirichlet parameter governing freq. T */
		GibbsParameterShPtr				kappaParam;				/**< updateable parameter representing kappa, the ts/tv rate ratio */
		GibbsParameterShPtr				rateVarParam;			/**< updateable parameter representing the variance parameter in a discrete gamma relative rates model (equal to inverse of the gamma shape parameter) */
		GibbsHyperParameterShPtr		edgelenHyperParam;		/**< updateable parameter representing the parameter governing the mean of the edge length prior */
		GibbsParameterShPtr				pInvarParam;			/**< updateable parameter representing the proportion of invariable sites parameter */
		GibbsParameterShPtr				meanRatesParam;			/**< updateable parameter representing the mean relative rate parameter */
		unsigned						ncat;					/**< number of discrete gamma rate categories */

		MCMCSettings					settings;				/**< pointer to settings structure holding options user chose for this run */
		unsigned						totalSteps;				/**< number of steps completed out of the MCMCSettings::niterations requested */
		unsigned						ntax;					/**< number of taxa in trees being explored */

		mutable NxsOutputStreamWrapperShPtr paramFPtr;					/**< file used for saving sampled parameter values */
#		if defined(DUMP_TREES)
			mutable std::ofstream			treef;					/**< file used for saving sampled trees */
#		endif
		Tree							*tree;					/**< tree object */
		double							prev_ln_likelihood;		/**< holds log likelihood from previous step */
		double							prev_ln_edgelen_prior;	/**< holds log edge length prior from previous step */
		double							prev_ln_polytomy_prior;	/**< holds log polytomy prior from previous step */
		double							prev_ln_prior;			/**< holds log prior from previous step */
		NxsOutputDestinationDescription splitSummaryDestination;
		PhoMCMC();
		
	public:
		void SetTaxaManager(PhoTaxaManager * inTaxaMgr)
			{
			taxaMgr = inTaxaMgr;
			}
		void SetTreesManager(PhoTreesManager * inTreesMgr)
			{
			treesMgr = inTreesMgr;
			}
		void SetCharactersManager(PhoCharactersManager * inCharactersMgr)
			{
			charactersMgr = inCharactersMgr;
			}
		PhoTaxaManager * taxaMgr;
		PhoTreesManager * treesMgr;
		PhoCharactersManager * charactersMgr;
#		if defined (CORBA_CLIENT_PHYCAS)
			boost::shared_ptr<cipresCORBA::PutTreeWrapper> treeConsumer;
#		endif

		~PhoMCMC();
		friend PhoMCMCCreator;
	};
#endif
