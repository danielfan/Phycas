#include "phycas/force_include.h"
#if defined(SUPPORT_GETTREES)

#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/cleartrees_settings.hpp"
#include "phycas/trees/trees_manager.hpp"
#include "ncl/command/nxs_set_cmd_param.hpp"
#include <boost/shared_ptr.hpp>
#include "ncl/command/nxs_command_manager.hpp"
#include <boost/bind.hpp>

void PhoTreesManager::SetupClearTrees(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<PhoTreesManager, ClearTreesSettings> PhycasClearTreesSettingsCmdCallback;
	typedef boost::shared_ptr<ClearTreesSettings> ClearTreesSettingsPtr;
	ClearTreesSettingsPtr settingsStruct = ClearTreesSettingsPtr(new ClearTreesSettings);
	PhycasClearTreesSettingsCmdCallback vClearTreesCmdCallback = PhycasClearTreesSettingsCmdCallback(this, &PhoTreesManager::HandleClearTrees, settingsStruct);

	typedef NxsAutoCommand<PhycasClearTreesSettingsCmdCallback> PhycasClearTreesCommand;
	typedef boost::shared_ptr<PhycasClearTreesCommand> PhycasClearTreesCommandPtr;

	PhycasClearTreesCommandPtr vClearTreesCommand = PhycasClearTreesCommandPtr(new PhycasClearTreesCommand("ClearTrees", vClearTreesCmdCallback));
	vClearTreesCommand->SetDescription("Removes the specified trees from memory");
	NxsIndexSetCmdOption * vtoClearCmdOpt = InstantiateTreeSetOption(std::string(), &settingsStruct->toClear, std::string(), this, true, kCmdPermBasicUser);
	vClearTreesCommand->AddUnnamedSetting(NxsCmdOptionShPtr(vtoClearCmdOpt), false);
	cmdMgr->AddCommand(vClearTreesCommand);
	}
#endif //if defined(SUPPORT_GETTREES)

