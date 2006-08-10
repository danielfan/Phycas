#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/taxa/nxs_taxa_manager.hpp"
#include "ncl/taxa/nxs_taxa_block.hpp"
#include "ncl/nxs_exception.hpp"
#include "ncl/misc/algorithm_extensions.hpp"
#include "ncl/misc/nxs_index_set_output.hpp"
#include "ncl/command/nxs_set_cmd_param.hpp"
#include "ncl/output/nxs_output.hpp"
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
	assert(0);
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
	