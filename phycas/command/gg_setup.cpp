#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/gg_settings.hpp"
#include "ncl/command/nxs_cmd_param.hpp"
#include <boost/shared_ptr.hpp>
#include "phycas/modules/gg/gg.hpp"
#include "ncl/command/nxs_command_manager.hpp"
#include <boost/bind.hpp>

void GG::SetupGG(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<GG, GGSettings> PhycasGGSettingsCmdCallback;
	typedef boost::shared_ptr<GGSettings> GGSettingsPtr;
	GGSettingsPtr settingsStruct = GGSettingsPtr(new GGSettings);
	PhycasGGSettingsCmdCallback vGGCmdCallback = PhycasGGSettingsCmdCallback(this, &GG::HandleGG, settingsStruct);

	typedef NxsAutoCommand<PhycasGGSettingsCmdCallback> PhycasGGCommand;
	typedef boost::shared_ptr<PhycasGGCommand> PhycasGGCommandPtr;

	PhycasGGCommandPtr vGGCommand = PhycasGGCommandPtr(new PhycasGGCommand("GG", vGGCmdCallback));
	vGGCommand->SetDescription("Computes Gelfand-Ghosh statistics for model selection (4-taxa only at this point)");
	NxsStringCmdOption * vTaxonsetCmdOpt = new NxsStringCmdOption("Taxonset", &settingsStruct->taxonset, "None", false, true, kCmdPermBasicUser);
	vGGCommand->AddKeyword(NxsCmdOptionShPtr(vTaxonsetCmdOpt), false);
	vTaxonsetCmdOpt->SetDescription("The name of the taxon set to be used for computing the Gelfand and Ghosh measures.");
	cmdMgr->AddCommand(vGGCommand);
	}
