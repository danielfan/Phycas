#include "phycas/force_include.h"
#if defined(READ_PHYCAS_BLOCK)
#include "phycas/floor.hpp"
#include "ncl/command/nxs_auto_command.hpp"
#include "ncl/misc/nxs_test_impl.hpp"
#include "phycas/rand/distribution_manager.hpp"
#include "phycas/modules/sampler/sampler.hpp"
#include "phycas/modules/mcmc/mcmc.hpp"
#include "phycas/misc/pwd.hpp"
#include "phycas/taxa/taxa_manager.hpp"
#include "phycas/trees/trees_manager.hpp"

#if defined(SCM_MODULE)//POL-23Nov2004
#	include "phycas/modules/scm/scm.hpp"
#   define FLOOR_SETUP_SCM_CMD	scm->SetupSCM(this)
#else
#   define FLOOR_SETUP_SCM_CMD
#endif

#if defined(GG_MODULE)//POL-14April2005
#	include "phycas/modules/gg/gg.hpp"
#   define FLOOR_SETUP_GG_CMD	gg->SetupGG(this)
#else
#   define FLOOR_SETUP_GG_CMD
#endif

#if defined(PHYLODIVERSITY_MODULE) //POL-3May04
#	include "phycas/modules/phylodiversity/pd.hpp"
#   define FLOOR_SETUP_PHYLODIVERSITY_CMD	Phylodiversity & pd = PhylodiversitySingletonHolder::Instance(); pd.SetTaxaManager(phoTaxaMgr.get()); pd.SetTreesManager(phoTreesMgr.get()); pd.SetupPhylodiversity(this);
#else
#   define FLOOR_SETUP_PHYLODIVERSITY_CMD
#endif

#if defined(SRQ_MODULE) //POL 9-September-2004
#	include "phycas/modules/srq/srq.hpp"
#   define FLOOR_SETUP_SRQ_CMD	srq->SetupSRQ(this)
#else
#   define FLOOR_SETUP_SRQ_CMD
#endif
#if defined(SUPPORT_GETTREES)
#   define FLOOR_SETUP_GETTREES_CMD SetupGetTrees(this)
#   define FLOOR_SETUP_CLEARTREES_CMD PhoTreesManager::GetInstance().SetupClearTrees(this)
#else
#   define FLOOR_SETUP_GETTREES_CMD
#   define FLOOR_SETUP_CLEARTREES_CMD
#endif

#if defined (WIN_PHOREST)
#   define FLOOR_SETUP_PWD_CMD PWDSingletonHolder::Instance().SetupPWD(this);
#else
#   define FLOOR_SETUP_PWD_CMD
#endif

void PhoFloor::AddRecognizedCommands()
	{
	SetupAliasCommand(this);
	SetupCopyingCommand(this);
	DistributionManagerSingletonHolder::Instance().SetupDistributionCommand(this);
	SetupExecuteCommand(this);
	SamplerSingletonHolder::Instance().SetupSampleCommand(this);
	//SetupHelpCommand(this);
	SetupQuitCommand(this);
	PhoMCMC & mcmc = MCMCSingletonHolder::Instance();
	mcmc.SetTaxaManager(phoTaxaMgr.get()); 
	mcmc.SetTreesManager(phoTreesMgr.get());
	mcmc.SetCharactersManager(phoCharactersMgr.get());
	mcmc.SetupMCMCCommand(this);
	FLOOR_SETUP_PWD_CMD;
	FLOOR_SETUP_SCM_CMD;
	FLOOR_SETUP_GG_CMD;
	FLOOR_SETUP_PHYLODIVERSITY_CMD;
	FLOOR_SETUP_SRQ_CMD;
	FLOOR_SETUP_CLEARTREES_CMD;
	FLOOR_SETUP_GETTREES_CMD;
	SetupSetCommand(this);
	//PhoTaxaManager & phoTaxaMgr = PhoTaxaManager::GetInstance();
	phoTaxaMgr->SetupTaxStatusCommand(this);
	phoTaxaMgr->SetupTaxSetCommand(this);
	phoTaxaMgr->SetupIncludeTaxaCommand(this);
	phoTaxaMgr->SetupExcludeTaxaCommand(this);
	SetAllowAbbreviations(true, true);
  	}
#endif // defined(READ_PHYCAS_BLOCK)
