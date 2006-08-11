#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/taxset_settings.hpp"
#include "ncl/command/nxs_command_manager.hpp"
#include "ncl/command/nxs_restricted_string_cmd_param.hpp"
#include "ncl/command/nxs_set_cmd_param.hpp"
#include "phycas/taxa/taxa_manager.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

void PhoTaxaManager::SetupTaxSetCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<PhoTaxaManager, TaxSetSettings> PhycasTaxSetSettingsCmdCallback;
	typedef boost::shared_ptr<TaxSetSettings> TaxSetSettingsPtr;
	TaxSetSettingsPtr settingsStruct = TaxSetSettingsPtr(new TaxSetSettings);
	PhycasTaxSetSettingsCmdCallback vTaxSetCmdCallback = PhycasTaxSetSettingsCmdCallback(this, &PhoTaxaManager::HandleTaxSet, settingsStruct);

	typedef NxsAutoCommand<PhycasTaxSetSettingsCmdCallback> PhycasTaxSetCommand;
	typedef boost::shared_ptr<PhycasTaxSetCommand> PhycasTaxSetCommandPtr;

	PhycasTaxSetCommandPtr vTaxSetCommand = PhycasTaxSetCommandPtr(new PhycasTaxSetCommand("TaxSet", vTaxSetCmdCallback));
	vTaxSetCommand->SetDescription("Specifies a name to be assigned to a group of taxa");
	RestrictNameCmdOption * vnameCmdOpt = new RestrictNameCmdOption(std::string(), &settingsStruct->name, std::string(), boost::bind(&PhoTaxaManager::GetReservedSetNames, this), true, kCmdPermBasicUser);
	NxsIndexSetCmdOption * vtaxSetCmdOpt = InstantiateTaxSetOption(std::string(), &settingsStruct->taxSet, std::string(), this, true, kCmdPermBasicUser);
	vTaxSetCommand->AddUnnamedSetting(NxsCmdOptionShPtr(vnameCmdOpt), false);
	vTaxSetCommand->ExpectEqualsSign();
	vTaxSetCommand->AddUnnamedSetting(NxsCmdOptionShPtr(vtaxSetCmdOpt), false);
	cmdMgr->AddCommand(vTaxSetCommand);
	}
