#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/excludetaxa_settings.hpp"
#include "phycas/taxa/taxa_manager.hpp"
#include "ncl/command/nxs_set_cmd_param.hpp"
#include <boost/shared_ptr.hpp>
#include "ncl/command/nxs_command_manager.hpp"
#include <boost/bind.hpp>

void PhoTaxaManager::SetupExcludeTaxaCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<PhoTaxaManager, ExcludeTaxaSettings> PhycasExcludeTaxaSettingsCmdCallback;
	typedef boost::shared_ptr<ExcludeTaxaSettings> ExcludeTaxaSettingsPtr;
	ExcludeTaxaSettingsPtr settingsStruct = ExcludeTaxaSettingsPtr(new ExcludeTaxaSettings);
	PhycasExcludeTaxaSettingsCmdCallback vExcludeTaxaCmdCallback = PhycasExcludeTaxaSettingsCmdCallback(this, &PhoTaxaManager::HandleExcludeTaxa, settingsStruct);

	typedef NxsAutoCommand<PhycasExcludeTaxaSettingsCmdCallback> PhycasExcludeTaxaCommand;
	typedef boost::shared_ptr<PhycasExcludeTaxaCommand> PhycasExcludeTaxaCommandPtr;

	PhycasExcludeTaxaCommandPtr vExcludeTaxaCommand = PhycasExcludeTaxaCommandPtr(new PhycasExcludeTaxaCommand("ExcludeTaxa", vExcludeTaxaCmdCallback));
	vExcludeTaxaCommand->SetDescription("Removes a set of taxa from further analyses, but does not alter thnumbering of taxa. The taxa can be reactivated by using the IncludeTaxa command");
	NxsIndexSetCmdOption * vtaxSetCmdOpt = InstantiateTaxSetOption(std::string(), &settingsStruct->taxSet, std::string(), this, true, kCmdPermBasicUser);
	vExcludeTaxaCommand->AddUnnamedSetting(NxsCmdOptionShPtr(vtaxSetCmdOpt), false);
	cmdMgr->AddCommand(vExcludeTaxaCommand);
	}
