#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/floor.hpp"

void PhoFloor::SetupCopyingCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsVoidExecuteCommandCallback<PhoFloor> PhycasCopyingSettingsCmdCallback;
	PhycasCopyingSettingsCmdCallback vCopyingCmdCallback = PhycasCopyingSettingsCmdCallback(this, &PhoFloor::HandleCopying);

	typedef NxsAutoCommand<PhycasCopyingSettingsCmdCallback> PhycasCopyingCommand;
	typedef boost::shared_ptr<PhycasCopyingCommand> PhycasCopyingCommandPtr;

	PhycasCopyingCommandPtr vCopyingCommand = PhycasCopyingCommandPtr(new PhycasCopyingCommand("Copying", vCopyingCmdCallback));
	vCopyingCommand->SetDescription("");
	cmdMgr->AddCommand(vCopyingCommand);
	}
