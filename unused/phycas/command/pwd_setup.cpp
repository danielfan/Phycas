#include "phycas/force_include.h"
#if defined(WIN_PHOREST)

#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/misc/pwd.hpp"
#include "ncl/command/nxs_command_manager.hpp"

void PWD::SetupPWD(NxsCommandManager *cmdMgr)
	{
	typedef NxsVoidExecuteCommandCallback<PWD> PhycasPWDSettingsCmdCallback;
	PhycasPWDSettingsCmdCallback vPWDCmdCallback = PhycasPWDSettingsCmdCallback(this, &PWD::HandlePWD);

	typedef NxsAutoCommand<PhycasPWDSettingsCmdCallback> PhycasPWDCommand;
	typedef boost::shared_ptr<PhycasPWDCommand> PhycasPWDCommandPtr;

	PhycasPWDCommandPtr vPWDCommand = PhycasPWDCommandPtr(new PhycasPWDCommand("PWD", vPWDCmdCallback));
	vPWDCommand->SetDescription("Displays the current working directory");
	cmdMgr->AddCommand(vPWDCommand);
	}
#endif //if defined(WIN_PHOREST)

