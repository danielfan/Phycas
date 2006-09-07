#ifndef PHO_TREE_MANAGER_H
#define PHO_TREE_MANAGER_H

#include "phypy/src/ncl/trees/nxs_trees_manager.hpp"
#include "phypy/src/oldphycas/tree_id.hpp"	// for the GetIDFromDescription function, perhaps we should move it.
class ClearTreesSettings;
class NxsCommandManager;
/*----------------------------------------------------------------------------------------------------------------------
|	NxsTreesManager with a programmatic interface for adding new tree descriptions and the ability to handle a 
|	TreeStatus command.
*/
class PhoTreesManager : 
  public NxsTreesManager
  	{
	public:
		CmdResult			AddNewTrees(const std::vector<FullTreeDescription> &v, unsigned d = UINT_MAX);
		CmdResult			HandleTreeStatus() const;
#		if defined(SUPPORT_GETTREES)
			void				SetupClearTrees(NxsCommandManager *);
			CmdResult			HandleClearTrees(ClearTreesSettings *);
			CmdResult			ClearTrees(const NxsIndexSet &toClear);
#		endif		
		STATIC_DATA_FUNC TreeID		GetIDFromDescription(const FullTreeDescription &);
	
		PhoTreesManager(PhoTaxaManager & taxaMgr);
		virtual ~PhoTreesManager(){}
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Calls NxsTreesManager::AddNewTrees and then displays the current tree status.
*/
inline CmdResult PhoTreesManager::AddNewTrees(const std::vector<FullTreeDescription> &v, unsigned d)
	{
	CmdResult r =  NxsTreesManager::AddNewTrees(v, d);
	if (r == kCmdSucceeded)
		r = HandleTreeStatus();
	return r;
	}
	
#endif

