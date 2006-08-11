#include "phycas/force_include.h"

// Macro below invalidates conditional likelihood arrays for entire tree before computing likelihood
// Do not turn on except to check likelihood calculations
//
//#define SAFE_MODE 

//The SAVE_FREQUENCY_TRIPLETS macro is for testing dirichlet prior on base frequencies
// When defined, four files are saved (ACG.txt, ACT.txt, AGT.txt and CGT.txt). In each
// of these files, tripets of are saved with each triplet comprising three normalized
// frequency parameters. For example, in ACG.txt, each triplet is
// (freqAParam/ACGtotal, freqCParam/ACGtotal, freqGParam/ACGtotal)
// where ACGtotal = freqAParam + freqCParam + freqGParam
// If a flat dirichlet prior is used, a barycentric plot of these triplets should 
// provide even coverage of the triangle. A 4-dimentional equivalent of a barycentric
// plot (i.e. pyramid rather than triangle) would be better.
//
//#define SAVE_FREQUENCY_TRIPLETS 

//#define USING_MARKS_PRIOR

// If OVERDISPERSE_INITIAL_VALUES is defined, one round of Gibbs sampling will be performed
// in PhoMCMC::StartingSampler() using the prior as the target distribution. This will choose 
// initial values for model parameters from their prior, causing each chain to start from
// values that are probably not all that good. Suggested if using Gelman-Rubin ratio to assess
// convergence.
//
//#define OVERDISPERSE_INITIAL_VALUES
#include "ncl/output/nxs_output.hpp"
#include "ncl/output/temp_basic_output_operators.hpp"
#include "ncl/misc/string_extensions.hpp"
#include "phycas/characters/characters_manager.hpp"
#include "phycas/misc/long_operation_manager.hpp"
#include "phycas/modules/mcmc/bush_master.hpp"
#include "phycas/modules/mcmc/edge_move.hpp"
#include "phycas/modules/mcmc/edge_len_prior.hpp" // to be deprecated
#include "phycas/modules/mcmc/gibbs_hyperparam.hpp"
#include "phycas/modules/mcmc/gibbs_param.hpp"
#include "phycas/modules/mcmc/larget_simon_move.hpp"
#include "phycas/modules/mcmc/mcmc.hpp"
#include "phycas/modules/mcmc/mcmc_settings_output.hpp"
#include "phycas/modules/mcmc/tree_scaler_move.hpp"
#include "phycas/rand/lot.hpp"
#include "phycas/taxa/taxa_manager.hpp"
#include "phycas/trees/draw_context.hpp"
#include "phycas/trees/hky_ad_hoc.hpp" // to be deprecated
#include "phycas/trees/split_manager.hpp"
#include "phycas/trees/topology_manager.hpp"
#include "phycas/trees/tree.hpp"
#include "phycas/trees/tree_output.hpp"
#if defined (CORBA_CLIENT_PHYCAS)
#	include "phycas/corba/put_tree_wrapper.hpp"
	using cipresCORBA::PutTreeWrapper;
#endif

using ncl::endl;
using ncl::flush;
using std::string;
using std::vector;

class MajRuleTreeSummary
	{
	public:
		double cutoffProportion;
		Tree tree;
		
		MajRuleTreeSummary(const double cutoffProp, const SplitManager &splitMgr, const UInt nTax)
			:cutoffProportion(cutoffProp)
			{
			splitMgr.BuildMajorityRuleTree(tree, cutoffProportion, nTax);
			splitMgr.SetBrlensFromSplits(tree);
			}
	};
	
template<class OUT_STREAM>
class GenericPrinterClass<kGraphicalTreeOutStyle, MajRuleTreeSummary, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const MajRuleTreeSummary & mrt)
			{
			outStream << mrt.cutoffProportion*100 << "% majority-rule consensus tree:\n";
			Emit<kGraphicalTreeOutStyle>(outStream, mrt.tree);
			}
	};
template<class OUT_STREAM>
class GenericPrinterClass<kNxsTreeTranslateCmdOutStyle, PhoTaxaManager, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const PhoTaxaManager & taxMgr)
			{
			const NxsIndexSet & nis = taxMgr.GetActiveSet();
			outStream << "translate\n" ;
			for (NxsIndexSet::const_iterator sIt = nis.begin(); sIt != nis.end(); ++sIt)
				{
				if (sIt != nis.begin())
					outStream << ",\n";
				outStream << *sIt << ' ' << GetAsNexusToken(taxMgr.GetLabel(*sIt)); 
				}
			outStream << ";\n";
			}
	};

template <class OUT_STREAM>
void PrintSplitSummaryLabels(OUT_STREAM & out);
template <>
void PrintSplitSummaryLabels<NxsTable>(NxsTable & table)
	{
	table.AddString("Split");
	table.SetLeftMargin();
	table.AddString("Length");
	table.AddString("Freq.");
	table.AddString("Prob.");
#	if defined(SHOW_SPLITID)
		table.AddString("Id");
#	endif
	table.SetRightMargin();
	table.HyphenLine();
	table.SetTopMargin();
	}
template <class OUT_STREAM>
void PrintSplitSummaryLabels(OUT_STREAM & out)
	{
	EmitFirstColumn(out, "Split") << NextCol << "Length" << NextCol << NextCol << "Freq." << NextCol << "Prob.";
#	if defined(SHOW_SPLITID)
		out << NextCol << "Id";
#	endif
	EndFirstRow(out);
	}
template<class OUT_STREAM>
void PrintSplitInfoTabular(OUT_STREAM & out, const SplitMap::const_iterator & it, const UInt totalTrees, const bool showTrivials, const double probCutoff, std::string &scratch);

template<>
void PrintSplitInfoTabular<NxsTable>(NxsTable & table, const SplitMap::const_iterator & it, const UInt totalTrees, const bool showTrivials, const double probCutoff, std::string & scratch)
	{	//@POL should fix isTrivial so that it is accurate even for immediate descendant of root node (which is a tip)
	const bool trivial_split = (it->second.IsTrivial() || it->first.CalcComplexity() == 1);
	const double prob = it->second.GetPosteriorProbability(totalTrees);
	if (prob >= probCutoff && (!trivial_split || showTrivials))
		{
		scratch.clear();
		it->first.CreateAndAppendPatternRepresentation(&scratch);
		table.AddString(scratch);
		table.AddDouble(it->second.GetMeanEdgeLength());
		table.AddUInt(it->second.GetFrequency());
		table.AddDouble(prob);
#		if defined(SHOW_SPLITID)
			table.AddString(it->first.CreateIdRepresentation());
#		endif
		}
	}

template<class OUT_STREAM>
void PrintSplitInfoTabular(OUT_STREAM & out, const SplitMap::const_iterator & i, const UInt totalTrees, const bool showTrivials, const double probCutoff, std::string &scratch)
	{ //@POL should fix isTrivial so that it is accurate even for immediate descendant of root node (which is a tip)
	const bool trivial_split = (i->second.IsTrivial() || i->first.CalcComplexity() == 1);
	const double prob = i->second.GetPosteriorProbability(totalTrees);
	if (prob >= probCutoff && (!trivial_split || showTrivials))
		{
		scratch.clear();
		i->first.CreateAndAppendPatternRepresentation(&scratch);
		EmitFirstColumn(out, scratch) << NextCol << i->second.GetMeanEdgeLength() << NextCol << i->second.GetFrequency() << NextCol << prob;
#		if defined(SHOW_SPLITID)
			out << NextCol << i->first.CreateIdRepresentation();
#		endif
		EndRow(out);
		}
	}
	
template<>
class GenericPrinterClass<kSplitSummaryTable, SplitManager, NxsWritableStream>
	{
	public:
		GenericPrinterClass(NxsWritableStream & outStream, const SplitManager & splitMgr)
		{
		try
			{
			NxsTable table;
			const UInt width = 12;// should be a parameter
			const UInt precision = 6;// should be a parameter
			NxsTableCell::def_width     = width;
			NxsTableCell::def_precision = precision;
			Emit<kSplitSummaryTable, SplitManager, NxsTable>(table, splitMgr);
			table.SetBottomMargin();
			table.Show(outStream);
			}
		catch(NxsTable::NxsX_InsufficientWidth & )
			{
			outStream << "  Table of splits too wide to show given current width";
			}
		}
	};

template<class OUT_STREAM>
class GenericPrinterClass<kSplitSummaryTable, SplitManager, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const SplitManager & splitMgr)
		{
		if (splitMgr.GetNumSplits() > 0)
			{
			const bool showTrivials = false; //@should be a parameter
			const double probCutoff = 0.0; //@should be a parameter
			std::string s;
			PrintSplitSummaryLabels<OUT_STREAM>(outStream);
			const SplitMap::const_iterator endIt = splitMgr.end();
			const UInt total = splitMgr.GetNumTreesRecorded();
			for (SplitMap::const_iterator i = splitMgr.begin(); i != endIt; ++i) 
				PrintSplitInfoTabular<OUT_STREAM>(outStream, i, total, showTrivials, probCutoff, s);
			}
		}
	};


/*----------------------------------------------------------------------------------------------------------------------
|	Constructs PhoMCMC object.
*/
PhoMCMC::PhoMCMC()
	:moveSchedule(MCMCSettings()),
	taxaMgr(NULL),
	treesMgr(NULL)
	{
	tree				= NULL;
	ntax				= 0;
	totalSteps			= 0;
	r					= LotShPtr(new Lot());
	splitMgr			= SplitManagerShPtr(new SplitManager(taxaMgr));
	topoMgr				= TopologyManagerShPtr(new TopologyManager(taxaMgr));
	ncat				= 1;
	}
		
/*----------------------------------------------------------------------------------------------------------------------
|	Destructs PhoMCMC object
*/
PhoMCMC::~PhoMCMC()
	{
	delete tree;
	}
		
/*----------------------------------------------------------------------------------------------------------------------
|	Handles MCMC command. Right now, this is little more than a transfer of what was in HandleCopying. Taking baby
|	steps...
*/
CmdResult PhoMCMC::HandleMCMC(
	MCMCSettings *s)	/* settings structure */
	{
	settings = *s;
	splitSummaryDestination = s->splitOut;
	try 
		{
		PrepareSampler();
		ShowInfo();
		RunSampler();
		}
	catch(...)
		{
		NxsOutput::GetNullCheckingStream<NxsOutput::kError>() << "\n(something bad happened and, unbelievably, Phycas ignored it)" << endl;
		return kCmdFailedSilent;
		}
	return kCmdSucceeded;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Repository for things that must be done before ShowInfo function is called.
*/
void PhoMCMC::PrepareSampler()
	{
	//assert(settings != NULL);

	if (settings.rnseed == 0)
		r->UseClockToSeedThenSpin();
	else
		r->SetSeed(settings.rnseed);

	if (tree == NULL)
		tree = new Tree();

	freqAParam			= GibbsParameterShPtr(new GibbsParameter(1.0));
	freqCParam			= GibbsParameterShPtr(new GibbsParameter(1.0));
	freqGParam			= GibbsParameterShPtr(new GibbsParameter(1.0));
	freqTParam			= GibbsParameterShPtr(new GibbsParameter(1.0));
	kappaParam			= GibbsParameterShPtr(new GibbsParameter(1.0));
	rateVarParam		= GibbsParameterShPtr(new GibbsParameter(2.0));
	//pInvarParam		= GibbsParameterShPtr(new GibbsParameter(0.0));
	//meanRatesParam	= GibbsParameterShPtr(new GibbsParameter(1.0));

	// Note: edgelenHyperParam cannot be assigned here (must wait until tree points to something)
	// but can go ahead and assign to edgelenHyperPrior

	// Specify prior distribution for the edge length mean
	// Inverse gamma facts:
	//   mean     = 1/[b(a-1)]
	//   variance = 1/[b^2 (a-1)^2 (a-2)]
	//   a = 2 + epsilon
	//   b = 1/[(mean)(1+epsilon)]
	//   variance = mean^2/epsilon
	// How to calculate a and b given mean and variance (e.g. mean=2, variance=9):
	//   epsilon = mean^2/variance      e.g. 2^2/9 = 4/9
	//   a       = 2 + epsilon          e.g. 2 + 4/9 = 22/9
	//   b       = 1/[mean*(1+epsilon)] e.g. 1/[2*(1+4/9)] = 9/26
	//
	//assert(settings.edgeLenHyperPriorVar > 0.0);
	double hyperprior_epsilon = pow(settings.edgeLenHyperPriorMean, 2.0)/settings.edgeLenHyperPriorVar;
	double hyperprior_a = 2.0 + hyperprior_epsilon;
	double hyperprior_tmp = settings.edgeLenHyperPriorMean*(1.0 + hyperprior_epsilon);
	//assert(hyperprior_tmp > 0.0);
	double hyperprior_b = 1.0/hyperprior_tmp;
	edgelenHyperPrior = ProbabilityDistributionShPtr(new InverseGammaDistribution(hyperprior_a, hyperprior_b));

	kappaParam->SetRandomNumberGenerator(r);
	kappaParam->SetUnitWidth(settings.initUnitWidth);
	kappaParam->SetValue(settings.kappaStartingValue);

	freqAParam->SetRandomNumberGenerator(r);
	freqAParam->SetUnitWidth(settings.initUnitWidth);
	freqAParam->SetValue(1.0);

	freqCParam->SetRandomNumberGenerator(r);
	freqCParam->SetUnitWidth(settings.initUnitWidth);
	freqCParam->SetValue(1.0);

	freqGParam->SetRandomNumberGenerator(r);
	freqGParam->SetUnitWidth(settings.initUnitWidth);
	freqGParam->SetValue(1.0);

	freqTParam->SetRandomNumberGenerator(r);
	freqTParam->SetUnitWidth(settings.initUnitWidth);
	freqTParam->SetValue(1.0);

	rateVarParam->SetRandomNumberGenerator(r);
	rateVarParam->SetUnitWidth(settings.initUnitWidth);
	rateVarParam->SetValue(1.0);

	// Specify prior distribution for base frequencies
	freqParamPrior = ProbabilityDistributionShPtr(new GammaDistribution(1.0, 1.0));
	freqAParam->SetPriorDistr(freqParamPrior);
	freqCParam->SetPriorDistr(freqParamPrior);
	freqGParam->SetPriorDistr(freqParamPrior);
	freqTParam->SetPriorDistr(freqParamPrior);

	// Specify prior distribution for kappa
	kappaPrior = ProbabilityDistributionShPtr(new ExponentialDistribution(1.0));
	kappaParam->SetPriorDistr(kappaPrior);

	// Specify prior distribution for the variance of gamma-distributed relative rates
	rateVarPrior = ProbabilityDistributionShPtr(new GammaDistribution(1.0, 1.0));
	rateVarParam->SetPriorDistr(rateVarPrior);

	edgeLenPriorCalc = EdgeLenPriorCalculatorShPtr(new EdgeLenPriorCalculator());
	edgelenPrior = ProbabilityDistributionShPtr(new ExponentialDistribution(0.1));
	edgeLenPriorCalc->SetProbabilityDistribution(edgelenPrior);

	edgelenHyperParam = GibbsHyperParameterShPtr(new GibbsHyperParameter(1.0, edgelenPrior, edgeLenPriorCalc, tree));
	edgelenHyperParam->SetPriorDistr(edgelenHyperPrior);
	edgelenHyperParam->SetRandomNumberGenerator(r);
	edgelenHyperParam->SetUnitWidth(settings.initUnitWidth);
	edgelenHyperParam->SetValue(1.0);

	hky = HKYAdHocEvaluatorShPtr(new HKYAdHocEvaluator(freqAParam, freqCParam, freqGParam, freqTParam, kappaParam, rateVarParam, pInvarParam, meanRatesParam, ncat));

	ncat = settings.ncat;
	hky->SetNumRateCategories(settings.ncat);
	hky->SetIgnoreData(settings.ignoreData);
	hky->SetTree(tree);
	hky->CopyData(*taxaMgr, *charactersMgr);
#   if defined(CORBA_CLIENT_PHYCAS)
		const VecString & taxLabels = hky->GetTaxonLabels();
		treeConsumer = boost::shared_ptr<PutTreeWrapper>(new PutTreeWrapper(taxLabels));
#   endif
	ntax = hky->GetNTaxa();

	Split::CalcStatics(ntax);
	try
		{
		TreeManip tm(tree);
		//tm.SimpleBuildStarTree(ntax, &r, 0.5);
		//tree->BuildTreeFromString("(a:0.3831192178245, b:0.3831192178245, (c:0.3831192178245, d:0.3831192178245):0.3831192178245)", true, false); // const char *s, bool translateNames, bool readAsRooted
		tm.SimpleBuildYuleTree(*r, 2.0, ntax);
		}
	catch(XBadTreeStructure & x)
		{
		NxsOutput::GetNullCheckingStream<NxsOutput::kError>()  << "XBadTreeStructure exception thrown: \n  " << x.msg << endl; //mthoutok
		throw x;
		}
	catch(XBadTreeDef & x)
		{
		NxsOutput::GetNullCheckingStream<NxsOutput::kError>()  << "XBadTreeDef exception thrown: \n  " << x.msg << endl; //mthoutok
		throw x;
		}

	moveSchedule.ResetMoves(tree, r);
	
	NxsOutputStream * outStream = NxsOutput::GetOutputStreamPtr();
	if (outStream != NULL)
		Emit<kGraphicalTreeOutStyle>(*outStream, *tree);
		// Create attribute objects for all nodes in the tree
	hky->PrepareTree();
	topoMgr->SetSplitMgr(splitMgr);

		// Set up the polytomy prior vector
		// Note: requires initialized BushMaster object
	topoPrior = PriorSummary(settings, ntax);
		// Note: vPolytomyPrior[M] now holds the log of the prior for the resolution class having M + 1 internal nodes
		//prev_ln_polytomy_prior = vPolytomyPrior[ntax - 3]; //@ bug? had been overwritten a few lines down by:
		//													//@ prev_ln_polytomy_prior = 0.0;
	
	prev_ln_likelihood = 0.0;
	if (!settings.ignoreData)
		{
#		if defined(SAFE_MODE)
			hky->InvalidateTree();	//@POL when everything is working, this line should be able to be deleted
#		endif
		prev_ln_likelihood = hky->CalcLogLikelihood();

#if 0
#	if defined(POL_UNDERFLOW_CORRECTION)
		std::ofstream doof("ufl_lnL.txt");
#	else
		std::ofstream doof("std_lnL.txt");
#	endif
		std::string tmp;
		StrPrintF(tmp, "%.8f", prev_ln_likelihood);
		doof << tmp << std::endl;
		doof << hky->GetNPatterns() << std::endl;
		hky->DebugShowPatterns(doof);
		doof.close();
		exit(0);
#endif
		}
	prev_ln_polytomy_prior = 0.0; // @ bug?  seems like it should be prev_ln_polytomy_prior = topoPrior.GetLnTopologyPrior(*tree); 
	prev_ln_edgelen_prior = edgeLenPriorCalc->GetLnEdgeLenPrior(*tree);
	prev_ln_prior = prev_ln_polytomy_prior + prev_ln_edgelen_prior;
	moveSchedule.Reset(settings);
	lastRecordedChainState.clear();
#   if defined(CHECK_WITH_PAUP)
		DebugCheckWithPAUP();
		exit(0);
#   endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a NEXUS file containing a TREES block with the starting tree and a PAUP block that computes the 
|	log-likelihood of the starting tree using the HKY85 model. The file contains a printed comment that provides the 
|	log-likelihood computed by this program for comparison.
*/
void PhoMCMC::DebugCheckWithPAUP()
	{
	assert(tree != NULL);
	string newick;
	tree->AppendNewickRepresentation(newick, true, true, true, 20);	// use 20 decimal places in edge lengths for accurate comparison

	// Open file for saving data, trees and a PAUP block
	//
	std::ofstream outf("debugcheckwithpaup.nex");
	outf.setf(std::ios::fixed, std::ios::floatfield); 
	outf.setf(std::ios::showpoint);
	outf << "#nexus\n" << std::endl;
	// Save the data in the form of a TAXA block/CHARACTERS block combination
	//
	hky->DebugSaveData(outf);

	// Save an initial PAUP block the sole purpose of which is to specify
	// that branch lengths encountered in TREES blocks are to be preserved
	//
	outf << "\nbegin paup;\n";
	outf << "  set storebrlens;\n";
	outf << "end;\n\n" << std::endl;
	// Save the starting tree in a TREES block
	//
	outf << "begin trees;" << std::endl;
	outf << "  translate" << std::endl;
	for (unsigned i = 0; i < ntax; ++i)
		{
		outf << "  " << i << " '" << taxaMgr->GetLabel(i) << '\'';
		if (i < ntax - 1)
			outf << ',' << std::endl;
		else
			outf << std::endl;
		}
	outf << "  ;" << std::endl;
	outf << "  utree startingtree = " << newick << ';' << std::endl;
	outf << "end;" << std::endl;

	outf << std::endl;

	// Save parameter values in the file in the form of a NEXUS printed comment
	//
	double tratio             = hky->CalcTRatio();
	double kappa              = hky->GetKappa();
	double rate_variance      = hky->GetRateVariance();
	double shape              = hky->GetGammaShape();
	std::vector<double> pi    = hky->GetNucleotideFrequencies();

	outf << "[!" << std::endl;
	outf << "Emprical base frequencies = (" << pi[0] << ',' << pi[1] << ',' << pi[2] << ',' << pi[3] << ')' << std::endl;
	outf << "Kappa                     = " << kappa << std::endl;
	outf << "Tratio                    = " << tratio << std::endl;
	outf << "Shape                     = " << shape << std::endl;
	outf << "RateVar                   = " << rate_variance << std::endl;
	outf << "Number of patterns        = " << hky->GetNPatterns() << std::endl;
	outf << "Log-likelihood            = " << prev_ln_likelihood << std::endl;
	outf << ']' << std::endl;

	outf << std::endl;

	// Save a PAUP block that calculates the likelihood under the HKY85 model for the 
	// starting tree
	//
	string s;
	outf << "begin paup;" << std::endl;
	outf << "  set criterion=likelihood;" << std::endl;

	if (ncat == 1)
		outf << StrPrintF(s, "  lset nst=2 basefreq=(%.8f %.8f %.8f) tratio=%.8f rates=equal;\n", pi[0], pi[1], pi[2], tratio);
	else
		outf << StrPrintF(s, "  lset nst=2 basefreq=(%.8f %.8f %.8f) tratio=%.8f rates=gamma ncat=%d shape=%.8f;\n", pi[0], pi[1], pi[2], tratio, ncat, shape);
	outf << "  lscores 1 / userbrlens;" << std::endl;
	outf << "end;" << std::endl;

#   if 0
		outf << '[' << std::endl;
		hky->DebugShowPatterns(outf);
		outf << ']' << std::endl;
#   endif

	outf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Informs user of what is about to happen, summarizing the MCMC run settings.
*/
void PhoMCMC::ShowInfo()
	{
	settings.rnseed = r->GetInitSeed();  //@ we might want to discuss whether we should try to propagate this change in seed settings back to ncl (we might not want to.  This only changes the settings if the seed was set to 0, and we might want to preserve the "use the clock" setting in this case)
	NxsOutputStream * outStream = NxsOutput::GetOutputStreamPtr();
	if (outStream != NULL)
		{
		*outStream << "\nMCMC command:";
		Emit<kVerboseOutStyle>(*outStream, settings);
		*outStream << "\n  Number of taxa:                     " << hky->GetNTaxa();
		*outStream << "\n  Number of sites:                    " << hky->GetNChar();
		*outStream << "\n  Number of patterns:                 " << hky->GetNPatterns();
		*outStream << endl;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls StartingSampler, then runs chain(s) for `MCMCSettings::niterations' steps, calling StoppingSampler when done.
|	SampleChainState is called every `sampleEvery' steps (assuming sampleEvery > 0), ReportChainState is called every 
|	`reportEvery' steps (assuming reportEvery > 0), and PlotBestTree is called every `plotEvery' steps (assuming 
|	plotEvery > 0).
*/
void PhoMCMC::RunSampler()
	{
	if (!StartingSampler())
		{
		// If here, means something went wrong (e.g. could not open output file)
		//@POL should there be a "run cancelled" message here?
		return;
		}
	const bool doSample	= (settings.sampleEvery > 0);
	const bool doReport	= (settings.reportEvery > 0);
	const bool doPlot   = (settings.plotEvery > 0);
	LongOperationManager & longOpMgr = LongOperationManagerSingletonHolder::Instance();
	for (totalSteps = 1L; totalSteps <= settings.niterations; ++totalSteps)
		{
		if (longOpMgr.CheckAbort(LongOperationManager::LongOp_MCMC))
			break;
		const int moveCategory = moveSchedule.GetMoveCategory(totalSteps);
		if (moveCategory & MoveSchedule::kMetropHastings)
			NextStep();
		if (moveCategory & MoveSchedule::kGibbsMove)
			GibbsSampling();
		const bool willSample = doSample && (totalSteps % settings.sampleEvery == 0);
		const bool willReport = doReport && (totalSteps % settings.reportEvery == 0);
		if (willSample || willReport)
			{
			const ChainState & cs = RefreshInternalCopyOfChainState();
			if (willSample)
				SampleChainState(cs);
			if (willReport)
				ReportChainState(cs);
			}
		if (doPlot && (totalSteps % settings.plotEvery == 0))
			PlotBestTree();
		}
	StoppingSampler();
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Propose next move in MCMC run.
*/
void PhoMCMC::NextStep()
	{
	assert(r);
	assert(tree);
	MCMCMove * currMoveAlias = moveSchedule.GetNextTreeMove(*r, *tree);
	assert(currMoveAlias != NULL);
	currMoveAlias->ProposeNewState();
	const double lnHastingsRatio = currMoveAlias->GetLnHastingsRatio();
	const double lnJacobian = currMoveAlias->GetLnJacobian();

	// Note: could compute prior for kappa here as follows:
	//
	//	double ln_kappa_prior = kappaDist->CalcKappaLnPrior(hky->GetKappa());
	//
	// however, this is not necessary. Because kappa does not change during
	// any of the possible moves, it will only cancel out and thus can be 
	// safely ignored
	//

	// Compute prior on edge lengths. First set mean of edgelenPrior to the current value
	// of the edge length mean parameter.
	//
	if (settings.useEdgeLenHyper)
		edgelenHyperParam->RefreshEdgeLenPrior();
	double ln_edgelen_prior  = edgeLenPriorCalc->GetLnEdgeLenPrior(*tree);

	double ln_polytomy_prior = 0.0;
#   if !defined(USING_MARKS_PRIOR) //@POL substitute topoPriorConstant for USING_MARKS_PRIOR when it becomes available
		if (!settings.topoPriorFlat)
			ln_polytomy_prior = topoPrior.GetLnTopologyPrior(*tree);
#   endif
	double ln_prior          = ln_edgelen_prior + ln_polytomy_prior;
	
	double ln_likelihood = 0.0;
	if (!settings.ignoreData)
		{
#		if defined(SAFE_MODE)
			hky->InvalidateTree();	//@POL when everything is working, this line should be able to be deleted
#		endif
		ln_likelihood = hky->CalcLogLikelihood();
		}

	double ln_accept_ratio = ln_likelihood + ln_prior - prev_ln_likelihood - prev_ln_prior + lnHastingsRatio + lnJacobian;

#   if defined(USING_MARKS_PRIOR)
		if (currMoveAlias == moveSchedule.bushMaster.get())		
			ln_accept_ratio += (moveSchedule.bushMaster->IsBirthMoveProposed() ? -1.0 : 1.0); //@POL substitute cTopoPrior for 1.0 when it is available
#   endif

	if (ln_accept_ratio >= 0.0 || std::log(r->Uniform()) <= ln_accept_ratio) //metrop-hastings-green acceptance rule
		{
		prev_ln_likelihood = ln_likelihood;
		prev_ln_edgelen_prior = ln_edgelen_prior;
		prev_ln_polytomy_prior = ln_polytomy_prior;
		prev_ln_prior = prev_ln_polytomy_prior + prev_ln_edgelen_prior;
		currMoveAlias->Accept();
		}
	else
		currMoveAlias->Revert();
	}

const ChainState & PhoMCMC::RefreshInternalCopyOfChainState()
	{   //@POL should ask model to write out values of all of its parameters
	lastRecordedChainState.totalSteps = totalSteps;
	lastRecordedChainState.lnLikelihood = prev_ln_likelihood;
	lastRecordedChainState.cBase[0] = freqAParam->GetValue();
	lastRecordedChainState.cBase[1] = freqCParam->GetValue();
	lastRecordedChainState.cBase[2] = freqGParam->GetValue();
	lastRecordedChainState.cBase[3] = freqTParam->GetValue();
	lastRecordedChainState.kappa = kappaParam->GetValue();
	lastRecordedChainState.rateVar = rateVarParam->GetValue();
	lastRecordedChainState.edgelenHyperPrior = edgelenHyperParam->GetValue();
	lastRecordedChainState.treeAlias = tree;
	return lastRecordedChainState;
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Tells the cold chain to sample its current state. It is called after every `settings.sampleEvery' calls of the
|	NextStep function.
*/
void PhoMCMC::SampleChainState(const ChainState & chainState) const
	{
	const unsigned m = tree->GetNNodes() - ntax;
	assert(m > 0);
	assert(m <= ntax - 2);
		// record the tree in memory
	splitMgr->RecordAllSplits(*chainState.treeAlias);
	topoMgr->AddTopology(*chainState.treeAlias);
	if (paramFPtr)
		EmitGenericPlotData(*paramFPtr, chainState) << endl;
#	if defined(DUMP_TREES)
		Emit<kNexusTreeCmdOutStyle>(treef, chainState) << std::endl;
#	endif
#   if defined (CORBA_CLIENT_PHYCAS)
		if (treeConsumer && chainState.treeAlias) 
			treeConsumer->PutTree(*chainState.treeAlias);
#   endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Uses slice sampling to obtain a sample from the full conditional distribution of one model parameter.
*/
void PhoMCMC::GibbsSampleParam(
  GibbsParameter *param,	/**< is the parameter to be sampled */
  bool ignoreData)			/**< indicates whether sampling from prior (true) or posterior (false) */
	{
	//@POL eventually want to replace hky argument with a functor, but problem is that we actually need two functors: 
	// one to compute the posterior density, the other to tell the model (i.e. hky in this case) when this particular
	// parameter's value has changed.
	//
	param->SetModel(hky);
	param->IgnoreData(ignoreData);
	param->GetNextSample(); 
	param->AdaptUnitWidth(settings.sliceDivisions); //@POL note using _kappa_ slice divisions here
	prev_ln_likelihood = param->GetLastSampledLnL();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Cycles through the substitution model parameters, updating each in turn using slice sampling.
*/
void PhoMCMC::GibbsSampling(
  bool ignore_data)	/**< is a flag that, if true, causes the slice sampler to explore the prior for each parameter */
	{
	// Hyperparameter needs to go first because it does not actually calculate the likelihood and thus 
	// prev_ln_likelihood will end up being reported as 0.0 if edgelenHyperParam is the last parameter
	// updated
	// 
	if (settings.useEdgeLenHyper)
		{
		GibbsSampleParam(edgelenHyperParam.get(), ignore_data);
		}

	GibbsSampleParam(freqAParam.get(), ignore_data);
	GibbsSampleParam(freqCParam.get(), ignore_data);
	GibbsSampleParam(freqGParam.get(), ignore_data);
	GibbsSampleParam(freqTParam.get(), ignore_data);
	GibbsSampleParam(kappaParam.get(), ignore_data);
	if (ncat > 1)
		GibbsSampleParam(rateVarParam.get(), ignore_data);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Displays basic information about an MCMC run in progress. It is called after every `settings.reportEvery' calls of
|	the NextStep function, displaying the number of steps, the number of topologies examined, the log posterior density
|	of the cold chain, and an estimate of the time to completion.
*/
void PhoMCMC::ReportChainState(const ChainState &) const
	{
	double pct_done = 0.1*((unsigned)(1000.0*(double)(totalSteps)/(double)settings.niterations));
	string tmp;
	if (settings.sampleEvery == 0)
		StrPrintF(tmp, "\n%6.1f%% done (burnin)", pct_done);
	else
		StrPrintF(tmp, "\n%6.1f%% done (sampling)", pct_done);
	SliceStats cAss = freqAParam->SummarizeDiagnostics();
	SliceStats cCss = freqCParam->SummarizeDiagnostics();
	SliceStats cGss = freqGParam->SummarizeDiagnostics();
	SliceStats cTss = freqTParam->SummarizeDiagnostics();
	SliceStats kappass = kappaParam->SummarizeDiagnostics();
	SliceStats ratevarss = rateVarParam->SummarizeDiagnostics();
	SliceStats edgelenmeanss = edgelenHyperParam->SummarizeDiagnostics();
	if (kappass.nsamples == 0)
		{
		StrPrintF(tmp, "\n        cA       =%g", freqAParam->GetValue());
		StrPrintF(tmp, "\n        cC       =%g", freqCParam->GetValue());
		StrPrintF(tmp, "\n        cG       =%g", freqGParam->GetValue());
		StrPrintF(tmp, "\n        cT       =%g", freqTParam->GetValue());
		StrPrintF(tmp, "\n        kappa    =%g", kappaParam->GetValue());
		StrPrintF(tmp, "\n        edgelen  =%g", edgelenHyperParam->GetValue());
		if (ncat > 1)
			StrPrintF(tmp, "\n        rateVar=%g", rateVarParam->GetValue());
		}
	else
		{
		StrPrintF(tmp, "\n        cA       =%g, width=%g, failed=%g, evals=%g, n=%d", cAss.value, cAss.width, cAss.failed, cAss.evals, freqAParam->GetNumSamples());
		StrPrintF(tmp, "\n        cC       =%g, width=%g, failed=%g, evals=%g, n=%d", cCss.value, cCss.width, cCss.failed, cCss.evals, freqCParam->GetNumSamples());
		StrPrintF(tmp, "\n        cG       =%g, width=%g, failed=%g, evals=%g, n=%d", cGss.value, cGss.width, cGss.failed, cGss.evals, freqGParam->GetNumSamples());
		StrPrintF(tmp, "\n        cT       =%g, width=%g, failed=%g, evals=%g, n=%d", cTss.value, cTss.width, cTss.failed, cTss.evals, freqTParam->GetNumSamples());
		StrPrintF(tmp, "\n        kappa    =%g, width=%g, failed=%g, evals=%g, n=%d", kappass.value, kappass.width, kappass.failed, kappass.evals, kappaParam->GetNumSamples());
		StrPrintF(tmp, "\n        edgelen  =%g, width=%g, failed=%g, evals=%g, n=%d", edgelenmeanss.value, edgelenmeanss.width, edgelenmeanss.failed, edgelenmeanss.evals, edgelenHyperParam->GetNumSamples());
		if (ncat > 1)
			StrPrintF(tmp, "\n        rateVar=%g, width=%g, failed=%g, evals=%g, n=%d", ratevarss.value, ratevarss.width, ratevarss.failed, ratevarss.evals, rateVarParam->GetNumSamples());
		}

#if defined(POL_UNDERFLOW_CORRECTION)
#	if defined(POL_TEMP)
		std::string s = hky->DebugShowBouncesVector();
		StrPrintF(tmp, "\n        vBounces=%s", s.c_str());
#	endif
#endif

	freqAParam->ResetDiagnostics();
	freqCParam->ResetDiagnostics();
	freqGParam->ResetDiagnostics();
	freqTParam->ResetDiagnostics();
	kappaParam->ResetDiagnostics();
	edgelenHyperParam->ResetDiagnostics();
	rateVarParam->ResetDiagnostics();
	NxsOutput::GetNullCheckingStream<NxsOutput::kOutput>() << tmp << endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Displays the majority rule consensus tree based on the current frequencies of splits. It is called after every 
|	`settings.plotEvery' calls of the NextStep function. Assumes splitMgr is non-NULL.
*/
void PhoMCMC::PlotBestTree()
	{
	assert(splitMgr != NULL);
	NxsOutputStream * outStream = NxsOutput::GetOutputStreamPtr();
	if (outStream != NULL)
		{
		MajRuleTreeSummary summary(0.5, *splitMgr, ntax);
		Emit<kGraphicalTreeOutStyle>(*outStream, summary) << ncl::endl;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Repository for things that must be done just before the MCMC run begins. This is where subcontractors are given a 
|	chance to open log files. If something goes wrong (e.g. output file cannot be opened), returns false; otherwise
|	returns true.
*/
bool PhoMCMC::StartingSampler()
	{

#   if defined(OVERDISPERSE_INITIAL_VALUES)
		// Do one round of Gibbs sampling to draw starting values of all parameters
		// from their respective prior distributions.
		//
		GibbsSampling(true);
#   endif

	paramFPtr = NxsOutput::GetGenericOutputStreamShPtr(settings.paramOut);
	if (paramFPtr)
		EmitGenericPlotLabels<ChainState *>(*paramFPtr, (ChainState *)NULL) << flush;
#	if defined(DUMP_TREES)
		//@POL need to do this the right way
		treef.open("sampled_trees.tre");
		treef << "#nexus\n\n";
		treef << "begin trees;\n\t";
		Emit<kNxsTreeTranslateCmdOutStyle>(treef, *taxaMgr) << std::endl;
#	endif
	LongOperationManagerSingletonHolder::Instance().StartLongOperation(LongOperationManager::LongOp_MCMC);
	NxsOutputStream * outStream = NxsOutput::GetOutputStreamPtr();
	if (outStream != NULL)
		{
		Emit<kVerboseModelParameterOutStyle>(*outStream, *hky);
		Emit(*outStream, topoPrior);
		Emit<kVerboseOutStyle>(*outStream, moveSchedule);
		*outStream << "\n\nstarting mcmc" << endl;
		}
	stopwatch.Start();
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Repository for things that must be done just after the MCMC run ends. Subcontractors are given the opportunity to 
|	close any log files they have opened.
*/
void PhoMCMC::StoppingSampler()
	{
	NxsOutputStream *outStream = NxsOutput::GetOutputStreamPtr();

	stopwatch.Stop();
#if defined(POL_NEW_STOPWATCH)
	if (outStream != NULL)
		*outStream << "\n<stopping>\n" << stopwatch.GetTimeUsed(5) << endl;
#else
	if (outStream != NULL)
		*outStream << "\n<stopping>\n" << stopwatch.ElapsedTime() << endl;
#endif
#	if defined(DUMP_TREES)
		treef << "end;" << std::endl;
		treef.close();
#	endif
	LongOperationManagerSingletonHolder::Instance().StopLongOperation(LongOperationManager::LongOp_MCMC);
	if (outStream != NULL)
		{
		*outStream << "\n\nParameter summary:" << endl;
		hky->ShowUpdatableParameters(*outStream);
		*outStream << endl;
		}
	//NxsOutputStreamWrapperShPtr splitOut = NxsOutput::GetGenericOutputStreamShPtr(splitSummaryDestination);
	//if (splitOut && splitMgr)
	//	{
	//	*splitOut << "Split summary:\n";
	//	EmitGeneric<kSplitSummaryTable>(*splitOut.get(), *splitMgr);
	//	}
	}
