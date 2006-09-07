//#include "phycas/force_include.h"
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/command/nxs_mixed_cmd_param.hpp"
#include <iterator>
#include "phypy/src/ncl/misc/algorithm_extensions.hpp"
using std::string;
#if !defined TEMPLATED_MIXED_CMD_OPTION

	VecString NxsMixedCmdOption::GetValidArgument()
		{
		VecString s;
		VecString c;
		for (VecNxsCmdOptionShPtr::iterator sovIt = subOptVector.begin(); sovIt != subOptVector.end(); ++sovIt)
			{
			if (sovIt != subOptVector.begin())
				s.push_back("|");
			c = (*sovIt)->GetValidArgument();
			copy(c.begin(), c.end(), back_inserter(s));
			}
		return s;
		}
		
	NxsMixedCmdOption::NxsMixedCmdOption(
	  const string &n, 
	  unsigned *optIndex,
	  unsigned defIndex,
	  VecNxsCmdOptionShPtr  vectorOfSubOption, 
	  bool persist, 
	  CmdPermissionLevel pLevel)
	  	:NxsCmdOption(n, false, persist, pLevel),
	  	subOptionIndex(optIndex),
	  	subOptVector(vectorOfSubOption),
	  	indexBefCmd(defIndex),
	  	defaultIndex(defIndex)
	  	{
	  	NXS_ASSERT(vectorOfSubOption.size() > 1);
	  	NXS_ASSERT(defIndex < vectorOfSubOption.size());
	  	*optIndex = defaultIndex;
	  	for (VecNxsCmdOptionShPtr::iterator sovIt = subOptVector.begin();  sovIt != subOptVector.end();  ++sovIt)
	  		{
	  		if ((*sovIt)->CheckIfOld())
	  			{
	  			ResetCheckIfOld(true);
	  			break;
	  			}
	  		}
	  	}

	string NxsMixedCmdOption::GetDisplayType(bool indef,
	  bool plural) const
		{
		VecNxsCmdOptionShPtr::const_iterator soIt = subOptVector.begin();
		string retVal = (*soIt++)->GetDisplayType(indef, plural);
		for (; soIt != subOptVector.end(); ++soIt)
			retVal << " or " << (*soIt)->GetDisplayType(indef, plural); 
		return retVal;
		}
			
	VecString NxsMixedCmdOption::GetLegalChoices() const 
		{
		VecString retVec;
		for (VecNxsCmdOptionShPtr::const_iterator soIt = subOptVector.begin(); soIt != subOptVector.end(); ++soIt)
			{
			VecString tempVec = (*soIt)->GetLegalChoices(); 
			if (tempVec.empty())
				retVec.push_back((*soIt)->GetDisplayType(true));
			else
				nxs_copy(tempVec.begin(), tempVec.end(), back_inserter(retVec));
			}
		return retVec;
		}	
		
	bool NxsMixedCmdOption::IsCurrentlyValid()
		{
		NXS_ASSERT(*subOptionIndex < subOptVector.size());
		return (subOptVector[*subOptionIndex]->IsCurrentlyValid() ? true : FlagError(subOptVector[*subOptionIndex]->GetErrorState(), subOptVector[*subOptionIndex]->GetErrorSnippet()));
		}
		
	void NxsMixedCmdOption::StorePreviousValue()
		{
		indexBefCmd = *subOptionIndex;
		std::for_each(subOptVector.begin(), subOptVector.end(), boost::mem_fn(&NxsCmdOption::PrepareToRead));
		}
		
	/*----------------------------------------------------------------------------------------------------------------------
	|	This function Tries to read the token stream with each of its sub options.
	|	If none can read the stream it _should_  set the errorstate for 
	|	the mixed cmd option to the least egregious error in the subcommands.
	|	Currently it just uses the last errorstate that isn't unrecognized or miss_wd.
	| 	If all sub options give an unrecognized or miss_wd error, unrecognized is used as the NxsMixedCmdOption error
	|	
	*/
	bool NxsMixedCmdOption::ReadValue(
	 NxsToken &token,	/* the stream of tokens that are being read */
	 bool equalsAlreadyRead) /* true if the equals sign has already been removed from the token stream (or is NOT expected) */
	  	{
	  	if (!equalsAlreadyRead && !EatEqualsThenAdvance(token))
	  		return false;
	  	*subOptionIndex = 0;
	  	NxsTokenizerState b = token.GetTokenizerState();
	  	NxsTokenizerState lowestErrorStreamPosition = b;
	  	unsigned lowestErrorOption = UINT_MAX;
	  	VecNxsCmdOptionShPtr::iterator soIt = subOptVector.begin();
		
	  	while (soIt != subOptVector.end())
	  		{
	  		if ((*soIt)->ReadValue(token, true))
	  			return WasValidRead();
	  		if ((*soIt)->GetErrorState() != NxsCmdOption::unrecognized && (*soIt)->GetErrorState() != NxsCmdOption::miss_wd)
				{
				lowestErrorOption = *subOptionIndex;
				lowestErrorStreamPosition = token.GetTokenizerState();
				}
	  		token.SeekTokenizerState(b);
	  		++soIt;
	  		(*subOptionIndex) += 1 ;
	  		}
	  	*subOptionIndex = 0;
	  	if (lowestErrorOption != UINT_MAX)
	  		{
	  		token.SeekTokenizerState(lowestErrorStreamPosition);
	  		FlagError(subOptVector[lowestErrorOption]->GetErrorState(), subOptVector[lowestErrorOption]->GetErrorSnippet());
	  		}
	  	else
	  		FlagError(NxsCmdOption::unrecognized);
	  	return false;
	  	}
	  	
	string NxsMixedCmdOption::GetCurrentValueAsString() const
		{
		NXS_ASSERT(*subOptionIndex < subOptVector.size());
		return subOptVector[*subOptionIndex]->GetCurrentValueAsString();
		}
			
	void NxsMixedCmdOption::RevertValueBecauseCommandFailed()
		{
		*subOptionIndex = indexBefCmd;
		std::for_each(subOptVector.begin(), subOptVector.end(), boost::mem_fn(&NxsCmdOption::RevertBecauseCommandFailed));
		}
		
	void NxsMixedCmdOption::ReturnValueToDefault()
		{
		*subOptionIndex = defaultIndex;
		std::for_each(subOptVector.begin(), subOptVector.end(), boost::mem_fn(&NxsCmdOption::ReturnToDefault));
		}


	bool NxsMixedCmdOption::AbbreviateChoices() 
		{
		bool ret = true;
		for (VecNxsCmdOptionShPtr::iterator soIt = subOptVector.begin(); soIt != subOptVector.end(); ++soIt)
			{
			if (!(*soIt)->AbbreviateChoices())
				ret = false;
			}
		return ret;
		}	
		

	  	
#endif

