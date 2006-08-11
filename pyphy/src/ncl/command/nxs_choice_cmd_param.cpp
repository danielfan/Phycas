//#include "phycas/force_include.h"
#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/command/nxs_choice_cmd_param.hpp"
using std::string;
/*----------------------------------------------------------------------------------------------------------------------
|	returns a | separated list of legal choices
*/
string NxsChoiceCmdOption::GetDisplayType(
  bool ,
  bool ) const
  	{
  	RefreshChoices();
  	string retStr;
  	for (VecString::iterator cIt = choices.begin(); cIt != choices.end(); ++cIt)
  		{
  		if (cIt != choices.begin())
  			retStr << '|';
  		retStr += *cIt;
  		}
  	return retStr;
  	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns false if allowAbbrevChoices is true, but unambiguous abbreviations weren't available.
|	Note that if Abbreviations aren't allowed, the choices are capitalized because IsValidChoice uses the 
|	IsCapAbbreviation
*/
bool NxsChoiceCmdOption::PrepareAbbreviations() const
	{
	if (allowAbbrevChoices)
		{
		VecString tempNames = choices;
		if (SetToShortestAbbreviation(tempNames))
			{
			choices = tempNames;
			return true;
			}
		}
	for_each(choices.begin(), choices.end(), &Capitalize);
	return (!allowAbbrevChoices);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	returns true if the s matches one of the legal choices.  If not, the appropriate errState flags are set.
|   called by CmdOption::WasValidRead()
*/
bool NxsChoiceCmdOption::IsValidChoice(
  const string &s)
  	{
	if (allowNumbers && IsAnUnsigned(s, choiceIndex))
		{
		if (*choiceIndex < choices.size())
			{ // @ using indexing starting at  here, users may prefer numbering to start at 1
			*currentValue = choices[*choiceIndex];
			return true;
			}
		}
	VecString::iterator choIt = find_if(choices.begin(), choices.end(), NStrCapAbbrevEquals(s));
	if (choIt != choices.end())
		{
		*currentValue = *choIt;	//if the user supplied an abbreviation, change it to the full word
		*choiceIndex = (unsigned) (choIt - choices.begin());
		return true;
		}
	*choiceIndex = UINT_MAX;
	if (choiceProvider)
		FlagError(NxsCmdOption::unrecog_labile);
	else
		FlagError(NxsCmdOption::unrecognized);
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|		
*/
NxsChoiceCmdOption::NxsChoiceCmdOption(
  const string &n, 
  unsigned *i,
  string *manipVal,
  string def, 
  const char *choiceStr,
  bool allowAbbreviations, 
  bool persist, 
  CmdPermissionLevel pLevel)
  	:NxsStringCmdOption(n, manipVal, def, true,  persist, pLevel),
  	choiceIndex(i),
  	allowAbbrevChoices(allowAbbreviations),
	allowNumbers(false)
  	{
  	choices = SplitString(choiceStr, '|');
  	*currentValue = def;
  	IsCurrentlyValid();
  	assert(errState == no_err);
  	}

/*----------------------------------------------------------------------------------------------------------------------
|		
*/
NxsChoiceCmdOption::NxsChoiceCmdOption(
  const string &n, 
  unsigned *i,
  string *manipVal,
  string def, 
  VecStringSource providerOfChoices,
  bool allowAbbreviations, 
  bool persist, 
  CmdPermissionLevel pLevel)
    :NxsStringCmdOption(n, manipVal, def, false, persist, pLevel),
  	choiceProvider(providerOfChoices),
  	choiceIndex(i),
	allowAbbrevChoices(allowAbbreviations),
	allowNumbers(false)
  	{
  	*currentValue = def;
  	IsCurrentlyValid();
  	}
	
/*----------------------------------------------------------------------------------------------------------------------
|		
*/
NxsChoiceCmdOption::NxsChoiceCmdOption(
  const string		  & n, 
  unsigned			  * i,
  string			  * manipVal,
  string				def, 
  const VecString     & choiceVector,
  bool					allowAbbreviations, 
  bool					persist, 
  CmdPermissionLevel	pLevel)
  	:NxsStringCmdOption(n, manipVal, def, false, persist, pLevel),
  	choiceIndex(i),
	choices(choiceVector),
  	allowAbbrevChoices(allowAbbreviations),
	allowNumbers(false)
  	{
  	*currentValue = def;
  	IsCurrentlyValid();
  	}
  		
