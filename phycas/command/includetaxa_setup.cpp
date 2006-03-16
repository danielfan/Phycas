#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/includetaxa_settings.hpp"
#include "phycas/taxa/taxa_manager.hpp"
#include "ncl/command/nxs_set_cmd_param.hpp"
#include <boost/shared_ptr.hpp>
#include "ncl/command/nxs_command_manager.hpp"
#include <boost/bind.hpp>

void PhoTaxaManager::SetupIncludeTaxaCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<PhoTaxaManager, IncludeTaxaSettings> PhycasIncludeTaxaSettingsCmdCallback;
	typedef boost::shared_ptr<IncludeTaxaSettings> IncludeTaxaSettingsPtr;
	IncludeTaxaSettingsPtr settingsStruct = IncludeTaxaSettingsPtr(new IncludeTaxaSettings);
	PhycasIncludeTaxaSettingsCmdCallback vIncludeTaxaCmdCallback = PhycasIncludeTaxaSettingsCmdCallback(this, &PhoTaxaManager::HandleIncludeTaxa, settingsStruct);

	typedef NxsAutoCommand<PhycasIncludeTaxaSettingsCmdCallback> PhycasIncludeTaxaCommand;
	typedef boost::shared_ptr<PhycasIncludeTaxaCommand> PhycasIncludeTaxaCommandPtr;

	PhycasIncludeTaxaCommandPtr vIncludeTaxaCommand = PhycasIncludeTaxaCommandPtr(new PhycasIncludeTaxaCommand("IncludeTaxa", vIncludeTaxaCmdCallback));
	vIncludeTaxaCommand->SetDescription("Includes taxa that are in memory, but have been made inactive by anExcludeTaxa command. The numbering of taxa is not affected by either IncludeTaxa or ExcludeTaxa.");
	NxsIndexSetCmdOption * vtaxSetCmdOpt = InstantiateTaxSetOption(std::string(), &settingsStruct->taxSet, std::string(), this, true, kCmdPermBasicUser);
	vIncludeTaxaCommand->AddUnnamedSetting(NxsCmdOptionShPtr(vtaxSetCmdOpt), false);
	cmdMgr->AddCommand(vIncludeTaxaCommand);
	}
