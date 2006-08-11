#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/mcmc_settings.hpp"
#include "ncl/command/nxs_primitive_cmd_param.hpp"
#include "ncl/command/nxs_command_manager.hpp"
#include "phycas/modules/mcmc/mcmc.hpp"
#include "ncl/command/nxs_file_cmd_param.hpp"
#include "ncl/misc/nxs_test_impl.hpp"
#include "phycas/taxa/taxa_manager.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

void PhoMCMC::SetupMCMCCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<PhoMCMC, MCMCSettings> PhycasMCMCSettingsCmdCallback;
	typedef boost::shared_ptr<MCMCSettings> MCMCSettingsPtr;
	MCMCSettingsPtr settingsStruct = MCMCSettingsPtr(new MCMCSettings);
	PhycasMCMCSettingsCmdCallback vMCMCCmdCallback = PhycasMCMCSettingsCmdCallback(this, &PhoMCMC::HandleMCMC, settingsStruct);

	typedef NxsAutoCommand<PhycasMCMCSettingsCmdCallback> PhycasMCMCCommand;
	typedef boost::shared_ptr<PhycasMCMCCommand> PhycasMCMCCommandPtr;

	PhycasMCMCCommandPtr vMCMCCommand = PhycasMCMCCommandPtr(new PhycasMCMCCommand("MCMC", vMCMCCmdCallback));
	vMCMCCommand->SetDescription("Starts an MCMC run using the data currently stored in memory");
	UIntCmdOption * vNIterationsCmdOpt = new UIntCmdOption("NIterations", &settingsStruct->niterations, 1100, 1, UINT_MAX, true, kCmdPermBasicUser);
	UIntCmdOption * vNChainsCmdOpt = new UIntCmdOption("NChains", &settingsStruct->nchains, 4, 1, UINT_MAX, true, kCmdPermBasicUser);
	UIntCmdOption * vSampleEveryCmdOpt = new UIntCmdOption("SampleEvery", &settingsStruct->sampleEvery, 10, 0, NamedNumberSource<unsigned>("NIterations", &settingsStruct->niterations), true, kCmdPermBasicUser);
	UIntCmdOption * vReportEveryCmdOpt = new UIntCmdOption("ReportEvery", &settingsStruct->reportEvery, 100, 0, NamedNumberSource<unsigned>("NIterations", &settingsStruct->niterations), true, kCmdPermBasicUser);
	UIntCmdOption * vPlotEveryCmdOpt = new UIntCmdOption("PlotEvery", &settingsStruct->plotEvery, 500, 0, NamedNumberSource<unsigned>("NIterations", &settingsStruct->niterations), true, kCmdPermBasicUser);
	NxsOutputCmdOption * vParamOutCmdOpt = new NxsOutputCmdOption("ParamOut", &settingsStruct->paramOut, NxsOutputDestinationDescription("samples.txt", false, true), true, kCmdPermBasicUser);
	NxsOutputCmdOption * vSplitOutCmdOpt = new NxsOutputCmdOption("SplitOut", &settingsStruct->splitOut, NxsOutputDestinationDescription(NxsOutputDestinationDescription::kRedirectToOutputStream), true, kCmdPermBasicUser);
	UIntCmdOption * vRNSeedCmdOpt = new UIntCmdOption("RNSeed", &settingsStruct->rnseed, 0, 0, UINT_MAX, true, kCmdPermBasicUser);
	BoolCmdOption * vIgnoreDataCmdOpt = new BoolCmdOption("IgnoreData", &settingsStruct->ignoreData, false, true, kCmdPermBasicUser);
	DblCmdOption * vBushMoveWeightCmdOpt = new DblCmdOption("BushMoveWeight", &settingsStruct->bushMoveWeight, 1, 0.0, DBL_MAX, true, kCmdPermBasicUser);
	UIntCmdOption * vGibbsEveryCmdOpt = new UIntCmdOption("GibbsEvery", &settingsStruct->gibbsEvery, 10, 0, NamedNumberSource<unsigned>("NIterations", &settingsStruct->niterations), true, kCmdPermBasicUser);
	UIntCmdOption * vSliceDivCmdOpt = new UIntCmdOption("SliceDiv", &settingsStruct->sliceDivisions, 3, 2, UINT_MAX, true, kCmdPermBasicUser);
	DblCmdOption * vInitUWidthCmdOpt = new DblCmdOption("InitUWidth", &settingsStruct->initUnitWidth, 50.0, DBL_MIN, DBL_MAX, true, kCmdPermBasicUser);
	DblCmdOption * vScalerMoveWeightCmdOpt = new DblCmdOption("ScalerMoveWeight", &settingsStruct->scalerMoveWeight, 0.0, 0.0, DBL_MAX, true, kCmdPermBasicUser);
	DblCmdOption * vLocalMoveWeightCmdOpt = new DblCmdOption("LocalMoveWeight", &settingsStruct->localMoveWeight, 2.0, 0.0, DBL_MAX, true, kCmdPermBasicUser);
	DblCmdOption * vKStartCmdOpt = new DblCmdOption("KStart", &settingsStruct->kappaStartingValue, 4.0, DBL_MIN, DBL_MAX, true, kCmdPermBasicUser);
	DblCmdOption * vKappaPriorMeanCmdOpt = new DblCmdOption("KappaPriorMean", &settingsStruct->kappaPriorMean, 1.0, 0, DBL_MAX, true, kCmdPermBasicUser);
	UIntCmdOption * vNCatCmdOpt = new UIntCmdOption("NCat", &settingsStruct->ncat, 1, 1, UINT_MAX, true, kCmdPermBasicUser);
	DblCmdOption * vPTopoPriorCmdOpt = new DblCmdOption("PTopoPrior", &settingsStruct->pTopoPrior, 0.5, 0.0, 1.0, true, kCmdPermBasicUser);
	BoolCmdOption * vFlatTopoPriorCmdOpt = new BoolCmdOption("FlatTopoPrior", &settingsStruct->topoPriorFlat, true, true, kCmdPermBasicUser);
	BoolCmdOption * vUseEdgeLenHyperPriorCmdOpt = new BoolCmdOption("UseEdgeLenHyperPrior", &settingsStruct->useEdgeLenHyper, true, true, kCmdPermBasicUser);
	DblCmdOption * vMeanEdgeLenHyperPriorCmdOpt = new DblCmdOption("MeanEdgeLenHyperPrior", &settingsStruct->edgeLenHyperPriorMean, 1, 0, DBL_MAX, true, kCmdPermBasicUser);
	DblCmdOption * vVarEdgeLenHyperPriorCmdOpt = new DblCmdOption("VarEdgeLenHyperPrior", &settingsStruct->edgeLenHyperPriorVar, 10, 0, DBL_MAX, true, kCmdPermBasicUser);
	BoolCmdOption * vResClassPriorCmdOpt = new BoolCmdOption("ResClassPrior", &settingsStruct->resClassPrior, false, true, kCmdPermBasicUser);
	DblCmdOption * vPolytomyLnPriorRatioCmdOpt = new DblCmdOption("PolytomyLnPriorRatio", &settingsStruct->polytomyLnPriorRatio, 0.0, -DBL_MAX, DBL_MAX, true, kCmdPermBasicUser);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vNIterationsCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vNChainsCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vSampleEveryCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vReportEveryCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vPlotEveryCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vParamOutCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vSplitOutCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vRNSeedCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vIgnoreDataCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vBushMoveWeightCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vGibbsEveryCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vSliceDivCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vInitUWidthCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vScalerMoveWeightCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vLocalMoveWeightCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vKStartCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vKappaPriorMeanCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vNCatCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vPTopoPriorCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vFlatTopoPriorCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vUseEdgeLenHyperPriorCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vMeanEdgeLenHyperPriorCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vVarEdgeLenHyperPriorCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vResClassPriorCmdOpt), false);
	vMCMCCommand->AddKeyword(NxsCmdOptionShPtr(vPolytomyLnPriorRatioCmdOpt), false);
	vNIterationsCmdOpt->SetDescription("Number of iterations to run the MCMC simulation");
	vNChainsCmdOpt->SetDescription("Number of chains to be run simultaneously (all but one of thesewill be heated chains)");
	vSampleEveryCmdOpt->SetDescription("Sample interval");
	vReportEveryCmdOpt->SetDescription("Report Interval");
	vPlotEveryCmdOpt->SetDescription("Frequency with which majority rule consensus tree will be displayedduring an MCMC run");
	vParamOutCmdOpt->SetDescription("Output destination for summaries of parameter values");
	vSplitOutCmdOpt->SetDescription("Destination for tabular summary of split frequencies");
	vRNSeedCmdOpt->SetDescription("Seed for pseudorandom number generator used in MCMC runs");
	vIgnoreDataCmdOpt->SetDescription("If specified, MCMC will run with no data and thus the posteriorwill be equivalent to the prior");
	vBushMoveWeightCmdOpt->SetDescription("Probability of proposing a dimension changing move to a lessresolved (or more resolved) tree topology");
	vGibbsEveryCmdOpt->SetDescription("Number of steps that must elapse before another round of Gibbssampling is done");
	vSliceDivCmdOpt->SetDescription("Number of divisions in slice after adapting slice sampler unit width.");
	vInitUWidthCmdOpt->SetDescription("Initial width of slice unit for sampling kappa values.");
	vScalerMoveWeightCmdOpt->SetDescription("Determines the probability of a tree scaling move, which scales alledges lengths in the tree up or down. Weights of all possible moves arenormalized to obtain probabilities.");
	vLocalMoveWeightCmdOpt->SetDescription("Determines probability of attempting a Larget-Simon LOCAL move.Weights of all possible moves are normalized to obtain the correspondingprobability of proposing a move of some type.");
	vKStartCmdOpt->SetDescription("Initial value of kappa parameter (should be moved to model description)");
	vKappaPriorMeanCmdOpt->SetDescription("Mean of the exponential prior for the kappa parameter (will makethis more general in the future)");
	vNCatCmdOpt->SetDescription("The number of rate categories for discrete gamma among site rate variation");
	vPTopoPriorCmdOpt->SetDescription("If TopoPriorFlat is false, PTopoPrior is used as the parameter ofthe geometric prior across resolution classes. A value near 0 means priorwill be close to flat across resolution classes (but not across trees), anda value near 1 produces a prior that strongly favors less-resolved trees.");
	vFlatTopoPriorCmdOpt->SetDescription("True specifies that prior should be flat across all topologies(regardless of the number of internal nodes). False specifies that priorshould be geometric, with PTopoPrior being the parameter.");
	vPolytomyLnPriorRatioCmdOpt->SetDescription("Log of the prior ratio for topology with m nodes to topology withm+1 nodes");
	StoredMsgSource sms_hasTaxa("there are no taxa in memory", std::string());
	NxsUIntCallbackWrapper left_hasTaxa(boost::bind(&PhoTaxaManager::GetNumTaxa, taxaMgr));
	NxsUIntValueWrapper right_hasTaxa(0);
	NxsTestShPtr test_hasTaxa = NxsTestShPtr(new UIntFuncToConstTest(left_hasTaxa, NxsTest::kGreaterThan, right_hasTaxa, sms_hasTaxa));
	vMCMCCommand->AddCommandAvailabilityReq(test_hasTaxa);
	cmdMgr->AddCommand(vMCMCCommand);
	}
