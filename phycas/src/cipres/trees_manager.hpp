/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#if defined(USING_OLDPHYCAS) || 1
#ifndef PHO_TREE_MANAGER_H
#define PHO_TREE_MANAGER_H

#include "phycas/src/ncl/trees/nxs_trees_manager.hpp"
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
			CmdResult			ClearTrees(const NxsIndexSet &toClear);
#		endif		
	
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
#endif //defined(USING_OLDPHYCAS)
