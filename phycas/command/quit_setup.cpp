#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/floor.hpp"

void PhoFloor::SetupQuitCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsVoidExecuteCommandCallback<PhoFloor> PhycasQuitSettingsCmdCallback;
	PhycasQuitSettingsCmdCallback vQuitCmdCallback = PhycasQuitSettingsCmdCallback(this, &PhoFloor::HandleQuit);

	typedef NxsAutoCommand<PhycasQuitSettingsCmdCallback> PhycasQuitCommand;
	typedef boost::shared_ptr<PhycasQuitCommand> PhycasQuitCommandPtr;

	PhycasQuitCommandPtr vQuitCommand = PhycasQuitCommandPtr(new PhycasQuitCommand("Quit", vQuitCmdCallback));
	vQuitCommand->SetDescription("");
	cmdMgr->AddCommand(vQuitCommand);
	}
