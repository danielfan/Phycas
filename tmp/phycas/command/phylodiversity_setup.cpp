#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/phylodiversity_settings.hpp"
#include "phycas/modules/phylodiversity/pd.hpp"
#include "ncl/command/nxs_command_manager.hpp"
#include "ncl/command/nxs_file_cmd_param.hpp"
#include "ncl/command/nxs_cmd_param.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

void Phylodiversity::SetupPhylodiversity(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<Phylodiversity, PhylodiversitySettings> PhycasPhylodiversitySettingsCmdCallback;
	typedef boost::shared_ptr<PhylodiversitySettings> PhylodiversitySettingsPtr;
	PhylodiversitySettingsPtr settingsStruct = PhylodiversitySettingsPtr(new PhylodiversitySettings);
	PhycasPhylodiversitySettingsCmdCallback vPhylodiversityCmdCallback = PhycasPhylodiversitySettingsCmdCallback(this, &Phylodiversity::HandlePhylodiversity, settingsStruct);

	typedef NxsAutoCommand<PhycasPhylodiversitySettingsCmdCallback> PhycasPhylodiversityCommand;
	typedef boost::shared_ptr<PhycasPhylodiversityCommand> PhycasPhylodiversityCommandPtr;

	PhycasPhylodiversityCommandPtr vPhylodiversityCommand = PhycasPhylodiversityCommandPtr(new PhycasPhylodiversityCommand("Phylodiversity", vPhylodiversityCmdCallback));
	vPhylodiversityCommand->SetDescription("Computes phylodiversity measures");
	NxsStringCmdOption * vTaxonsetCmdOpt = new NxsStringCmdOption("Taxonset", &settingsStruct->taxonset, "None", false, true, kCmdPermBasicUser);
	NxsOutputCmdOption * vPDOutputCmdOpt = new NxsOutputCmdOption("PDOutput", &settingsStruct->pdOut, NxsOutputDestinationDescription(), true, kCmdPermBasicUser);
	vPhylodiversityCommand->AddKeyword(NxsCmdOptionShPtr(vTaxonsetCmdOpt), false);
	vPhylodiversityCommand->AddKeyword(NxsCmdOptionShPtr(vPDOutputCmdOpt), false);
	vTaxonsetCmdOpt->SetDescription("The name of the taxon set to be treated as the focal group forpurposes of computing phylodiversity measures.");
	vPDOutputCmdOpt->SetDescription("Output destination for a tab-delimited summary of the outputcompatible with plotting with GnuPlot or importing into Excel.");
	cmdMgr->AddCommand(vPhylodiversityCommand);
	}
