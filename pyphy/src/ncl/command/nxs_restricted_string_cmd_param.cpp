//#include "phycas/force_include.h"
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/command/nxs_restricted_string_cmd_param.hpp"
using std::string;
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
NxsWordCmdOption::NxsWordCmdOption(
  const string &n, 
  string *manipVal,
  string def, 
  bool 		validityIsLabile,
  bool 		persist, 
  CmdPermissionLevel pLevel)
  	:NxsStringCmdOption(n, manipVal, def, validityIsLabile, persist, pLevel)
  	{
  	}
	
/*----------------------------------------------------------------------------------------------------------------------
|		
*/
RestrictNameCmdOption::RestrictNameCmdOption(
  const string &n, 
  string *manipVal,
  string def, 
  VecStringSource tabooProvider,
  bool persist, 
  CmdPermissionLevel pLevel)
  	:NxsWordCmdOption(n, manipVal, def, true, persist, pLevel),
  	tabooListProvider(tabooProvider)
  	{
  	NXS_ASSERT(tabooListProvider != NULL);
  	RefreshTabooList();
  	}
	
/*----------------------------------------------------------------------------------------------------------------------
|		
*/
RestrictNameCmdOption::RestrictNameCmdOption(
  const string & n, 
  string	  * manipVal,
  string		def, 
  VecString 	tabooVector,
  bool			persist, 
  CmdPermissionLevel pLevel)
  	:NxsWordCmdOption(n, manipVal, def, false, persist, pLevel),
  	tabooList(tabooVector)
  	{
  	*manipVal = def;
  	}	

/*----------------------------------------------------------------------------------------------------------------------
|		
*/
RestrictNameCmdOption::RestrictNameCmdOption(
  const string & n, 
  string	  * manipVal,
  string		def, 
  const char  * tabooStr,
  bool			persist, 
  CmdPermissionLevel pLevel)
  	:NxsWordCmdOption(n, manipVal, def, false, persist, pLevel)
  	{
  	tabooList = SplitString(tabooStr, '|');
  	*manipVal = def;
  	}
  	
/*----------------------------------------------------------------------------------------------------------------------
|	If the setting is a valid word, this function verifies that it is not in the list of restricted words
*/
bool RestrictNameCmdOption::IsCurrentlyValid()
	{
	if (!NxsWordCmdOption::IsCurrentlyValid())
		return false;
	RefreshTabooList();
  	VecString::iterator tabooWordIt = find_if(tabooList.begin(), tabooList.end(), NStrCaseInsensitiveEquals(*currentValue));
	return (tabooWordIt == tabooList.end() ? true : FlagError(NxsCmdOption::reserved_word));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	verifies that the current value of the setting is not a punctuation character (note: only flags single character
|	strings as errors.  this is not a test of whether or not the value contains a punctuation character.
*/
bool NxsWordCmdOption::IsCurrentlyValid()
	{
	return (IsLegalNexusWord(*currentValue) ? true :  FlagError(NxsCmdOption::illegal_punctuation));
	}

								
							
