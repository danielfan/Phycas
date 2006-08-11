#include "phycas/force_include.h"
#if defined(SUPPORT_GETTREES)

#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "phycas/command/gettrees_settings.hpp"
#include "phycas/trees/trees_manager.hpp"
#include "phycas/floor.hpp"
#include "ncl/command/nxs_set_cmd_param.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

void PhoFloor::SetupGetTrees(NxsCommandManager *cmdMgr)
	{
	typedef NxsExecuteCommandCallback<PhoFloor, GetTreesSettings> PhycasGetTreesSettingsCmdCallback;
	typedef boost::shared_ptr<GetTreesSettings> GetTreesSettingsPtr;
	GetTreesSettingsPtr settingsStruct = GetTreesSettingsPtr(new GetTreesSettings);
	PhycasGetTreesSettingsCmdCallback vGetTreesCmdCallback = PhycasGetTreesSettingsCmdCallback(this, &PhoFloor::HandleGetTrees, settingsStruct);

	typedef NxsAutoCommand<PhycasGetTreesSettingsCmdCallback> PhycasGetTreesCommand;
	typedef boost::shared_ptr<PhycasGetTreesCommand> PhycasGetTreesCommandPtr;

	PhycasGetTreesCommandPtr vGetTreesCommand = PhycasGetTreesCommandPtr(new PhycasGetTreesCommand("GetTrees", vGetTreesCmdCallback));
	vGetTreesCommand->SetDescription("Reads NEXUS Trees blocks in a file");
	NxsIndexSetCmdOption * vToSaveCmdOpt = InstantiateTreeSetOption("ToSave", &settingsStruct->toSave, std::string(), &PhoTreesManager::GetInstance(), true, kCmdPermBasicUser);
	vGetTreesCommand->AddKeyword(NxsCmdOptionShPtr(vToSaveCmdOpt), false);
	cmdMgr->AddCommand(vGetTreesCommand);
	}
#endif //if defined(SUPPORT_GETTREES)

