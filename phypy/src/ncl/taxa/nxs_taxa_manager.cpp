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
#include "phypy/src/ncl/taxa/nxs_taxa_block.hpp"
#include "phypy/src/ncl/nxs_exception.hpp"
#include "phypy/src/ncl/misc/algorithm_extensions.hpp"
#include "phypy/src/ncl/misc/nxs_index_set_output.hpp"
#include "phypy/src/ncl/command/nxs_set_cmd_param.hpp"
#include "phypy/src/ncl/output/nxs_output.hpp"
using std::string;
NxsTaxaManager::NxsTaxaManager()
	:reservedTaxAndSetNames(SplitString("ALL|INCLUDED|ACTIVE|REMAINDER", '|')),
	warnBeforeClearingTaxa(true)
	{
	}

NxsIndexSetCmdOption *InstantiateTaxSetOption(const string &n, NxsIndexSet *manip, const string &d, NxsTaxaManager *tmPtr, bool  persist, CmdPermissionLevel pLevel)
	{
	return InstantiateSetOptionGeneric<NxsTaxaManager>(n, manip, d, tmPtr, persist , "Tax", pLevel);
	}

CmdResult NxsTaxaManager::NewBlockRead(NxsBlock *newBlock)
	{
	NxsTaxaBlock * newTaxa = dynamic_cast<NxsTaxaBlock *>(newBlock);
	if (newTaxa != NULL)
		return ReplaceTaxa(newTaxa->GetTaxa());
	NXS_ASSERT(0);
	return kCmdSucceeded;
	}

CmdResult NxsTaxaManager::ReplaceTaxa(const VecString &v)
	{
	if (!taxLabels.empty())
		{
		if (warnBeforeClearingTaxa)
			{
			if (!NxsOutput::GetUserQueryPtr() || !NxsOutput::GetUserQueryPtr()->AskUserYesNoQuery("Replace Taxa", "Reading in a new taxa block will clear old taxa, trees and characters.  Do you want to proceed?"))
				return kCmdFailedSilent;
			}
		Clear();
		}
	AppendTaxa(v);
	return kCmdSucceeded;
	}
	
void NxsTaxaManager::DisplayTaxSet(const NxsIndexSet &s, bool compact)
	{
	NxsOutputStream *outS = NxsOutput::GetOutputStreamPtr();
	if (outS != NULL)
		{
		NxsIndexSetOutput so(s, *this, kTaxaSetOutput);
		if (compact)
			Emit<kConciseSetOutStyle>(*outS, so);
		else
			Emit<kVerboseOutStyle>(*outS, so);
		}
	}

void NxsTaxaManager::DisplayTaxSets(bool compact)
	{
	const VecString setNames = GetSetNames();
	for (VecString::const_iterator snIt = setNames.begin(); snIt != setNames.end(); ++snIt)
		DisplayTaxSet(*GetSet(*snIt), compact);
	}

NxsTaxaBlockAlias NxsTaxaManager::GetTaxaBlockReader()
	{
	SharedTaxaBlockPtr temp = SharedTaxaBlockPtr(new NxsTaxaBlock(*this));
	taxaBlocks.push_back(temp);
	return NxsTaxaBlockAlias(temp);	
	}	
	
