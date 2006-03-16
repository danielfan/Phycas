#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/execute_settings.hpp"
#include "phycas/floor.hpp"
#include "ncl/command/nxs_file_cmd_param.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

void PhoFloor::SetupExecuteCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<PhoFloor, ExecuteSettings> PhycasExecuteSettingsCmdCallback;
	typedef boost::shared_ptr<ExecuteSettings> ExecuteSettingsPtr;
	ExecuteSettingsPtr settingsStruct = ExecuteSettingsPtr(new ExecuteSettings);
	PhycasExecuteSettingsCmdCallback vExecuteCmdCallback = PhycasExecuteSettingsCmdCallback(this, &PhoFloor::HandleExecute, settingsStruct);

	typedef NxsAutoCommand<PhycasExecuteSettingsCmdCallback> PhycasExecuteCommand;
	typedef boost::shared_ptr<PhycasExecuteCommand> PhycasExecuteCommandPtr;

	PhycasExecuteCommandPtr vExecuteCommand = PhycasExecuteCommandPtr(new PhycasExecuteCommand("Execute", vExecuteCmdCallback));
	vExecuteCommand->SetDescription("Reads the specified NEXUS file, executing the commands from all recognized blocks");
	NxsInFileCmdOption * vFilenameCmdOpt = new NxsInFileCmdOption("Filename", &settingsStruct->fileToExecute, NxsInFilePath(std::string()), true, kCmdPermBasicUser);
	vExecuteCommand->AddKeyword(NxsCmdOptionShPtr(vFilenameCmdOpt), false);
	cmdMgr->AddCommand(vExecuteCommand);
	}
