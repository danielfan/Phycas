//#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/taxa/nxs_taxa_manager.hpp"
#include "ncl/trees/nxs_trees_manager.hpp"
#include "ncl/trees/nxs_tree_listener.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|	Stores alias to the NxsTreeListener and calls the NxsTreeListener::AddListener
*/
NxsTreeListener::NxsTreeListener(NxsTreesManager *t)
	:treeMgr(t)
	{
	if (treeMgr != NULL)
		treeMgr->AddListener(this);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	removes itself from the old NxsTreeListener list of listeners (via NxsTreeListener::RemoveListener) and add itself to 
|	the new one	Manager's list (NxsTreeListener::AddListener)
*/
void NxsTreeListener::ChangeManager(NxsTreesManager *t)
	{
	if (treeMgr != NULL)
		treeMgr->RemoveListener(this);
	treeMgr = t;
	if (treeMgr != NULL)
		treeMgr->AddListener(this);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	removes itself from the NxsTreesManager list of listeners (via NxsTreesManager::RemoveListener)
*/
NxsTreeListener::~NxsTreeListener()
	{
	if (treeMgr != NULL)
		treeMgr->RemoveListener(this);
	}
