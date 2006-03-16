#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/set_settings.hpp"
#include "phycas/floor.hpp"
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include "ncl/command/nxs_primitive_cmd_param.hpp"

void PhoFloor::SetupSetCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<PhoFloor, SetSettings> PhycasSetSettingsCmdCallback;
	typedef boost::shared_ptr<SetSettings> SetSettingsPtr;
	SetSettingsPtr settingsStruct = SetSettingsPtr(new SetSettings);
	PhycasSetSettingsCmdCallback vSetCmdCallback = PhycasSetSettingsCmdCallback(this, &PhoFloor::HandleSet, settingsStruct);

	typedef NxsAutoCommand<PhycasSetSettingsCmdCallback> PhycasSetCommand;
	typedef boost::shared_ptr<PhycasSetCommand> PhycasSetCommandPtr;

	PhycasSetCommandPtr vSetCommand = PhycasSetCommandPtr(new PhycasSetCommand("Set", vSetCmdCallback));
	vSetCommand->SetDescription("Sets general program behavior or changes settings of other commands withoutexecuting the command");
	UIntCmdOption * vOutWidthCmdOpt = new UIntCmdOption("OutWidth", &settingsStruct->outputWidth, 80, 1, UINT_MAX, true, kCmdPermBasicUser);
	vSetCommand->AddKeyword(NxsCmdOptionShPtr(vOutWidthCmdOpt), false);
	vOutWidthCmdOpt->SetDescription("Number of characters per line in the output display");
	cmdMgr->AddCommand(vSetCommand);
	}
