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

//#include "phycas/force_include.h"
#include "phycas/src/ncl/trees/nxs_trees_manager.hpp"
#include "phycas/src/ncl/taxa/nxs_taxa_manager.hpp"
#include "phycas/src/ncl/trees/nxs_trees_block.hpp"
#include <boost/bind.hpp>
#include "phycas/src/ncl/command/nxs_set_cmd_param.hpp"
#include "phycas/src/ncl/output/nxs_output.hpp"
using std::vector;
using std::string;
#define MTH_MARCH_3_2005_DEBUG

FullTreeDescription   NxsTreesManager::GetTree(unsigned i) const
	{
	if (i >= treeDescriptions.size())
		throw NxsX_IndexNotFound();
	return treeDescriptions[i];
	}
	
VecOfTreeDescriptions  NxsTreesManager::GetTrees(const NxsIndexSet &s) const
	{
	VecOfTreeDescriptions treeDescs;
	for (NxsIndexSet::const_iterator sIt = s.begin(); sIt != s.end(); ++sIt)
		{
		if (*sIt >= treeDescriptions.size())
			break;
		treeDescs.push_back(treeDescriptions[*sIt]);
		}
	return treeDescs;
	}
		
NxsIndexSetCmdOption *InstantiateTreeSetOption(const string &n, NxsIndexSet *manip, const string &d, NxsTreesManager *tmPtr, bool  persist, CmdPermissionLevel pLevel)
	{
	return InstantiateSetOptionGeneric<NxsTreesManager>(n, manip, d, tmPtr, persist , "Tree", pLevel);
	}


CmdResult NxsTreesManager::NewBlockRead(NxsBlock *newBlock)
	{
	NxsTreesBlock * newTrees = dynamic_cast<NxsTreesBlock *>(newBlock); //tucson2005: VC7 warning C4541: 'dynamic_cast' used on polymorphic type 'NxsBlock' with /GR-; unpredictable behavior may result nxs_taxa_listener.cpp
	if (newTrees != NULL)
		return AddNewTrees(newTrees->GetTreeDescriptions(), newTrees->GetDefaultTree());
	NXS_ASSERT(0);
	return kCmdSucceeded;
	}

void NxsTreesManager::AlertListeners(NxsTreeListener::TreeChangeType ct)
	{
	for_each(listeners.begin(), listeners.end(), boost::bind(&NxsTreeListener::TreesChanged, _1, this, ct));
	}

NxsTreesManager::NxsTreesManager(NxsTaxaManager & inTaxaMgr)
	:NxsTaxaListener(dynamic_cast<BaseTaxaManager*>(&inTaxaMgr)),
  	defaultTree(UINT_MAX),
	taxaMgr(inTaxaMgr)
	{
	};

NxsTreesManager::~NxsTreesManager()
	{
	}
	
CmdResult NxsTreesManager::AddNewTrees(const vector<FullTreeDescription> &v, unsigned def)
	{
	if (v.empty())
		return kCmdSucceeded;
	if (def != UINT_MAX)
		defaultTree = (unsigned)treeDescriptions.size() + def;
	treeDescriptions.reserve(treeDescriptions.size() + v.size());
	bool somePartialBrLens = false;
//	Split::CalcStatics(taxaMgr.GetSize()); //@shouldn't be necessary

	for (vector<FullTreeDescription>::const_iterator vIt = v.begin(); vIt != v.end(); ++vIt)
		{
		if (vIt->hasEdgeLens == FullTreeDescription::kSomeEdgeLengths)
			{
			FullTreeDescription tempD = *vIt;
			tempD.RemoveEdgeLengths();
			treeDescriptions.push_back(tempD);
			somePartialBrLens = true;
			}
		else
			treeDescriptions.push_back(*vIt);
		
#		if 0 //defined(MTH_MARCH_3_2005_DEBUG)
			Tree t;
			t.BuildTreeFromDescription(*vIt);
			t.DebugShowPreorder(std::cerr);
			for (UInt i = 0; i < 4; ++i)
				{
				std::cerr << "\n\n\nRerooting at " << i << "\n";
				t.RerootAtTip(i); //POL 10-Mar-2005 Mark, RerootAtTip now has a refresh_splits parameter that is true by default
				t.DebugShowPreorder(std::cerr);
				}
#		endif
		}
	if (somePartialBrLens)
		NxsOutput::GetNullCheckingStream<NxsOutput::kOutput>() << "Note:  some branch lengths were missing in the tree description.  No branch lengths have been stored for these trees.\n";	
	IndicesIncreased(GetSize(), true);
	AlertListeners(NxsTreeListener::kTreesAdded);
	return kCmdSucceeded;
	}


TreesBlockShPtr NxsTreesManager::GetTreesBlockReader()
	{
	TreesBlockShPtr temp = TreesBlockShPtr(new NxsTreesBlock(taxaMgr, *this)); //@@ shouldn't use the getInstance!!
	treeBlocks.push_back(temp);
	return temp;	
	}	
	

