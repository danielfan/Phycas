/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

//#include "phycas/force_include.h"
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/taxa/nxs_taxa_manager.hpp"
#include "phypy/src/ncl/trees/nxs_trees_manager.hpp"
#include "phypy/src/ncl/trees/nxs_tree_listener.hpp"

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
