#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/sample_settings.hpp"
#include "ncl/command/nxs_primitive_cmd_param.hpp"
#include "ncl/command/nxs_command_manager.hpp"
#include "phycas/modules/sampler/sampler.hpp"
#include "phycas/rand/distribution_manager.hpp"
#include "ncl/command/nxs_file_cmd_param.hpp"
#include "phycas/rand/distribution_command_param.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

void PhoSampler::SetupSampleCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<PhoSampler, SampleSettings> PhycasSampleSettingsCmdCallback;
	typedef boost::shared_ptr<SampleSettings> SampleSettingsPtr;
	SampleSettingsPtr settingsStruct = SampleSettingsPtr(new SampleSettings);
	PhycasSampleSettingsCmdCallback vSampleCmdCallback = PhycasSampleSettingsCmdCallback(this, &PhoSampler::HandleSample, settingsStruct);

	typedef NxsAutoCommand<PhycasSampleSettingsCmdCallback> PhycasSampleCommand;
	typedef boost::shared_ptr<PhycasSampleCommand> PhycasSampleCommandPtr;

	PhycasSampleCommandPtr vSampleCommand = PhycasSampleCommandPtr(new PhycasSampleCommand("Sample", vSampleCmdCallback));
	vSampleCommand->SetDescription("Samples random deviates from the specified probability distribution");
	DistributionCmdOption * vDistributionCmdOpt = new DistributionCmdOption("Distribution", &settingsStruct->distribution, "Uniform(0.0,1.0)", DistributionCmdOption::kDiscreteOrContinuous, 0, boost::bind(&DistributionManager::GetDistribution, &DistributionManagerSingletonHolder::Instance(), _1), true, kCmdPermBasicUser);
	UIntCmdOption * vNSampledCmdOpt = new UIntCmdOption("NSampled", &settingsStruct->nsampled, 100, 1, UINT_MAX, true, kCmdPermBasicUser);
	BoolCmdOption * vShowCmdOpt = new BoolCmdOption("Show", &settingsStruct->show, false, true, kCmdPermBasicUser);
	UIntCmdOption * vSeedCmdOpt = new UIntCmdOption("Seed", &settingsStruct->seed, 0, 0, UINT_MAX, true, kCmdPermBasicUser);
	NxsOutputCmdOption * vOutputDestinationCmdOpt = new NxsOutputCmdOption("OutputDestination", &settingsStruct->outDestination, NxsOutputDestinationDescription(), true, kCmdPermBasicUser);
	UIntCmdOption * vPrecisionCmdOpt = new UIntCmdOption("Precision", &settingsStruct->precision, 6, 0, 30, true, kCmdPermBasicUser);
	vSampleCommand->AddKeyword(NxsCmdOptionShPtr(vDistributionCmdOpt), false);
	vSampleCommand->AddKeyword(NxsCmdOptionShPtr(vNSampledCmdOpt), false);
	vSampleCommand->AddKeyword(NxsCmdOptionShPtr(vShowCmdOpt), false);
	vSampleCommand->AddKeyword(NxsCmdOptionShPtr(vSeedCmdOpt), false);
	vSampleCommand->AddKeyword(NxsCmdOptionShPtr(vOutputDestinationCmdOpt), false);
	vSampleCommand->AddKeyword(NxsCmdOptionShPtr(vPrecisionCmdOpt), false);
	vDistributionCmdOpt->SetDescription("Allows specification of the probability distribution from whichrandom deviates are to be sampled (either the distribution description orthe name of a previously specified distribution can be provided).");
	vNSampledCmdOpt->SetDescription("Number of random deviates to sample");
	vShowCmdOpt->SetDescription("if true, shows all random deviates sampled followed by a summary(e.g. mean and standard deviation); only the summary is output if NOSHOW is specified");
	vSeedCmdOpt->SetDescription("Random number generator seed (0 to use the system clock)");
	vOutputDestinationCmdOpt->SetDescription("Destination for output of sampled deviates");
	vPrecisionCmdOpt->SetDescription("Specifies precision (number of decimal places) to be used whendisplaying sampled values");
	cmdMgr->AddCommand(vSampleCommand);
	}
