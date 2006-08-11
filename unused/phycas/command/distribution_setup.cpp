#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/distribution_settings.hpp"
#include "phycas/rand/distribution_manager.hpp"
#include "ncl/command/nxs_command_manager.hpp"
#include "phycas/rand/distribution_command_param.hpp"
#include "ncl/command/nxs_restricted_string_cmd_param.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

void DistributionManager::SetupDistributionCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<DistributionManager, DistributionSettings> PhycasDistributionSettingsCmdCallback;
	typedef boost::shared_ptr<DistributionSettings> DistributionSettingsPtr;
	DistributionSettingsPtr settingsStruct = DistributionSettingsPtr(new DistributionSettings);
	PhycasDistributionSettingsCmdCallback vDistributionCmdCallback = PhycasDistributionSettingsCmdCallback(this, &DistributionManager::HandleNewDistribution, settingsStruct);

	typedef NxsAutoCommand<PhycasDistributionSettingsCmdCallback> PhycasDistributionCommand;
	typedef boost::shared_ptr<PhycasDistributionCommand> PhycasDistributionCommandPtr;

	PhycasDistributionCommandPtr vDistributionCommand = PhycasDistributionCommandPtr(new PhycasDistributionCommand("Distribution", vDistributionCmdCallback));
	vDistributionCommand->SetDescription("Defines a probability distribution. The name can then be used in anycommand that requires a distribution (for example when priors are specified).");
	RestrictNameCmdOption * vnameCmdOpt = new RestrictNameCmdOption(std::string(), &settingsStruct->name, std::string(), boost::bind(&DistributionManager::GetReservedNames, this), true, kCmdPermBasicUser);
	DistributionCmdOption * vdistribCmdOpt = new DistributionCmdOption(std::string(), &settingsStruct->distrib, std::string(), DistributionCmdOption::kDiscreteOrContinuous, 0, boost::bind(&DistributionManager::GetDistribution, this, _1), true, kCmdPermBasicUser);
	vDistributionCommand->AddUnnamedSetting(NxsCmdOptionShPtr(vnameCmdOpt), false);
	vDistributionCommand->ExpectEqualsSign();
	vDistributionCommand->AddUnnamedSetting(NxsCmdOptionShPtr(vdistribCmdOpt), false);
	vnameCmdOpt->SetDescription("Name to be assigned to the distribution");
	vdistribCmdOpt->SetDescription("Definition of the distribution");
	cmdMgr->AddCommand(vDistributionCommand);
	}
