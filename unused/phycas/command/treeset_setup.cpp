#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/treeset_settings.hpp"
#include "ncl/command/nxs_command_manager.hpp"
#include "ncl/command/nxs_restricted_string_cmd_param.hpp"
#include "ncl/command/nxs_set_cmd_param.hpp"
#include "phycas/trees/trees_manager.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

void PhoTreesManager::SetupTreeSetCommand(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<PhoTreesManager, TreeSetSettings> PhycasTreeSetSettingsCmdCallback;
	typedef boost::shared_ptr<TreeSetSettings> TreeSetSettingsPtr;
	TreeSetSettingsPtr settingsStruct = TreeSetSettingsPtr(new TreeSetSettings);
	PhycasTreeSetSettingsCmdCallback vTreeSetCmdCallback = PhycasTreeSetSettingsCmdCallback(this, &PhoTreesManager::HandleTreeSet, settingsStruct);

	typedef NxsAutoCommand<PhycasTreeSetSettingsCmdCallback> PhycasTreeSetCommand;
	typedef boost::shared_ptr<PhycasTreeSetCommand> PhycasTreeSetCommandPtr;

	PhycasTreeSetCommandPtr vTreeSetCommand = PhycasTreeSetCommandPtr(new PhycasTreeSetCommand("TreeSet", vTreeSetCmdCallback));
	vTreeSetCommand->SetDescription("Specifies a name to be assigned to a group of trees");
	RestrictNameCmdOption * vnameCmdOpt = new RestrictNameCmdOption(std::string(), &settingsStruct->name, std::string(), boost::bind(&PhoTreesManager::GetReservedSetNames, this), true, kCmdPermBasicUser);
	NxsIndexSetCmdOption * vtreeSetCmdOpt = InstantiateTreeSetOption(std::string(), &settingsStruct->treeSet, std::string(), this, true, kCmdPermBasicUser);
	vTreeSetCommand->AddUnnamedSetting(NxsCmdOptionShPtr(vnameCmdOpt), false);
	vTreeSetCommand->ExpectEqualsSign();
	vTreeSetCommand->AddUnnamedSetting(NxsCmdOptionShPtr(vtreeSetCmdOpt), false);
	cmdMgr->AddCommand(vTreeSetCommand);
	}
