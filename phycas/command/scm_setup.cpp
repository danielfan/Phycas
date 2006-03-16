#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/scm_settings.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include "phycas/modules/scm/scm.hpp"
#include "ncl/command/nxs_command_manager.hpp"
#include "ncl/command/nxs_file_cmd_param.hpp"

void SCM::SetupSCM(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<SCM, SCMSettings> PhycasSCMSettingsCmdCallback;
	typedef boost::shared_ptr<SCMSettings> SCMSettingsPtr;
	SCMSettingsPtr settingsStruct = SCMSettingsPtr(new SCMSettings);
	PhycasSCMSettingsCmdCallback vSCMCmdCallback = PhycasSCMSettingsCmdCallback(this, &SCM::HandleSCM, settingsStruct);

	typedef NxsAutoCommand<PhycasSCMSettingsCmdCallback> PhycasSCMCommand;
	typedef boost::shared_ptr<PhycasSCMCommand> PhycasSCMCommandPtr;

	PhycasSCMCommandPtr vSCMCommand = PhycasSCMCommandPtr(new PhycasSCMCommand("SCM", vSCMCmdCallback));
	vSCMCommand->SetDescription("Computes strict consensus merger supertree for trees currently stored");
	NxsOutputCmdOption * vTreeOuputCmdOpt = new NxsOutputCmdOption("TreeOuput", &settingsStruct->scmfile, NxsOutputDestinationDescription("scm.tre", false, false), true, kCmdPermBasicUser);
	vSCMCommand->AddKeyword(NxsCmdOptionShPtr(vTreeOuputCmdOpt), false);
	vTreeOuputCmdOpt->SetDescription("Output destination for strict consensus merger supertree description");
	cmdMgr->AddCommand(vSCMCommand);
	}
