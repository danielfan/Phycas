#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/taxstatus_settings.hpp"
#include "phycas/taxa/taxa_manager.hpp"
#include "ncl/command/nxs_primitive_cmd_param.hpp"
#include <boost/shared_ptr.hpp>
#include "ncl/command/nxs_command_manager.hpp"
#include <boost/bind.hpp>

void PhoTaxaManager::SetupTaxStatusCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<PhoTaxaManager, TaxStatusSettings> PhycasTaxStatusSettingsCmdCallback;
	typedef boost::shared_ptr<TaxStatusSettings> TaxStatusSettingsPtr;
	TaxStatusSettingsPtr settingsStruct = TaxStatusSettingsPtr(new TaxStatusSettings);
	PhycasTaxStatusSettingsCmdCallback vTaxStatusCmdCallback = PhycasTaxStatusSettingsCmdCallback(this, &PhoTaxaManager::HandleTaxStatus, settingsStruct);

	typedef NxsAutoCommand<PhycasTaxStatusSettingsCmdCallback> PhycasTaxStatusCommand;
	typedef boost::shared_ptr<PhycasTaxStatusCommand> PhycasTaxStatusCommandPtr;

	PhycasTaxStatusCommandPtr vTaxStatusCommand = PhycasTaxStatusCommandPtr(new PhycasTaxStatusCommand("TaxStatus", vTaxStatusCmdCallback));
	vTaxStatusCommand->SetDescription("Displays information about the current taxa");
	BoolCmdOption * vFullDisplayCmdOpt = new BoolCmdOption("FullDisplay", &settingsStruct->fullDisplay, true, true, kCmdPermBasicUser);
	BoolCmdOption * vShowExcludedCmdOpt = new BoolCmdOption("ShowExcluded", &settingsStruct->showExcluded, true, true, kCmdPermBasicUser);
	vTaxStatusCommand->AddKeyword(NxsCmdOptionShPtr(vFullDisplayCmdOpt), false);
	vTaxStatusCommand->AddKeyword(NxsCmdOptionShPtr(vShowExcludedCmdOpt), false);
	vFullDisplayCmdOpt->SetDescription("If true, a table with taxon number, name and inclusion status isshown. If false, a brief summary is displayed.");
	vShowExcludedCmdOpt->SetDescription("Controls whether or not excluded are listed in the summary");
	cmdMgr->AddCommand(vTaxStatusCommand);
	}
