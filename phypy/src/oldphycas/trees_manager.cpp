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
#include "phypy/src/oldphycas/taxa_manager.hpp"
#include "phypy/src/oldphycas/trees_manager.hpp"
#include "phypy/src/oldphycas/tree.hpp"
#include "phypy/src/ncl/output/nxs_output.hpp"
#include "phypy/src/oldphycas/cleartrees_settings.hpp"
using ncl::flush;

#undef CONSTRUCT_TREES_ON_INPUT

PhoTreesManager::PhoTreesManager(PhoTaxaManager & taxaMgr)
	:NxsTreesManager(taxaMgr)
	{} //wrap_scm.cpp is the only non-singelton instantiation

#if defined(SUPPORT_GETTREES)
	CmdResult PhoTreesManager::HandleClearTrees(ClearTreesSettings *s)
		{
		NXS_ASSERT(s != NULL);
		return ClearTrees(s->toClear);
		}
	CmdResult PhoTreesManager::ClearTrees(const NxsIndexSet & s)
		{
		return kCmdSucceeded;
		}
#endif


TreeID PhoTreesManager::GetIDFromDescription(const FullTreeDescription &ftd)
	{
	Tree refTree;
	refTree.BuildTreeFromDescription(ftd);
	refTree.RefreshID();
	return refTree.GetID();	
	}

/*----------------------------------------------------------------------------------------------------------------------
|	displays info on currently known trees to the user.
*/
CmdResult PhoTreesManager::HandleTreeStatus() const
	{
	if (NxsOutput::GetOutputStreamPtr() != NULL)
		{
		const unsigned ntrees = GetNumTrees();
		NxsOutputStream & outStream = *NxsOutput::GetOutputStreamPtr();
		if (ntrees == 1)
			outStream << "There is currently 1 tree.\n";
		else
			outStream << "There are currently " << ntrees << " trees.\n";
#		if defined(CONSTRUCT_TREES_ON_INPUT)		
			for (unsigned i = 0; i < ntrees; ++i)
				{
				outStream << '\t' << (i+1) << '\t' << treeDescriptions[i].BriefReport() << '\n' ;
				Tree t;
				try	{
					t.BuildTreeFromDescription(treeDescriptions[i]);
					//outStream << t;
					// debugging  t.ForAllLeaves(boost::bind(&Tree::DebugRerootingTest, &t, _1, boost::ref(outStream)));
					}
				catch  (XBadTreeDef & x)
					{
					outStream << x.msg  << '\n';
					}
				catch  (XBadTreeStructure & x)
					{
					outStream << x.msg  << '\n';
					}
				}
#		endif
		outStream << flush;
		}
	return kCmdSucceeded;
	}

#if defined(EFRON_BOOT_CMD_WRITER)
	/*----------------------------------------------------------------------------------------------------------------------
	|	
	*/
	CmdResult PhoTreesManager::ClearTrees(const NxsIndexSet &toClear)
		{
		if (toClear.empty())
			return kCmdSucceeded;
		NxsIndexSet::const_iterator indIt = toClear.end();		--indIt;
		bool someTreesWereRemoved = false;
		for (;;)
			{
			if (*indIt < GetNumTrees())
				{
				someTreesWereRemoved = true;
				treeDescriptions.erase(treeDescriptions.begin() + *indIt);
				if (defaultTree != UINT_MAX)
					{
					if (defaultTree > *indIt)
						--defaultTree;
					else if (defaultTree == *indIt)
						defaultTree = UINT_MAX;
					}
				}
			if (indIt != toClear.begin())
				--indIt;
			else
				break;
			}
		if (someTreesWereRemoved)
			{
			IndicesDecreased(GetSize());
			AlertListeners(NxsTreeListener::kTreesDeleted);
			}
		return kCmdSucceeded;
		}
#endif
