#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/alias_settings.hpp"
#include "phycas/floor.hpp"
#include "ncl/command/nxs_cmd_param.hpp"
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include "ncl/command/nxs_restricted_string_cmd_param.hpp"

void PhoFloor::SetupAliasCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<PhoFloor, AliasSettings> PhycasAliasSettingsCmdCallback;
	typedef boost::shared_ptr<AliasSettings> AliasSettingsPtr;
	AliasSettingsPtr settingsStruct = AliasSettingsPtr(new AliasSettings);
	PhycasAliasSettingsCmdCallback vAliasCmdCallback = PhycasAliasSettingsCmdCallback(this, &PhoFloor::HandleAlias, settingsStruct);

	typedef NxsAutoCommand<PhycasAliasSettingsCmdCallback> PhycasAliasCommand;
	typedef boost::shared_ptr<PhycasAliasCommand> PhycasAliasCommandPtr;

	PhycasAliasCommandPtr vAliasCommand = PhycasAliasCommandPtr(new PhycasAliasCommand("Alias", vAliasCmdCallback));
	vAliasCommand->SetDescription("Specifies a shortened name for a command or series of commands (note theexpanded command must be in single quotes)");
	RestrictNameCmdOption * valiasKeyCmdOpt = new RestrictNameCmdOption(std::string(), &settingsStruct->aliasKey, std::string(), boost::bind(&PhoFloor::GetCommandNames, this), true, kCmdPermBasicUser);
	NxsStringCmdOption * vexpansionCmdOpt = new NxsStringCmdOption(std::string(), &settingsStruct->expansion, std::string(), false, true, kCmdPermBasicUser);
	vAliasCommand->AddUnnamedSetting(NxsCmdOptionShPtr(valiasKeyCmdOpt), false);
	vAliasCommand->ExpectEqualsSign();
	vAliasCommand->AddUnnamedSetting(NxsCmdOptionShPtr(vexpansionCmdOpt), false);
	cmdMgr->AddCommand(vAliasCommand);
	}
