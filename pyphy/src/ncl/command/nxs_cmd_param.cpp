//#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_cmd_param.hpp"
#include "ncl/command/nxs_command.hpp"
#include "ncl/misc/utilities.hpp"
using std::string;

bool NxsCmdOption::IsAvailable() const
	{
	if (!HasPermission(NxsCommand::userPermLevel))
		return false;
	return PerformAllTests(availableRequirements, NULL);
	}

std::string NxsCmdOption::GetReasonForUnavailability() const
	{
	if (!HasPermission(NxsCommand::userPermLevel))
		return UserLevelTooLowString(GetPermissionLevel());
	std::string errM;
	PerformAllTests(availableRequirements, & errM);
	return errM;
	}

NxsRequiredToken::NxsRequiredToken(
  const string &rt)
  	:NxsCmdOption(string(), false, false, kCmdPermBasicUser),
  	reqToken(rt)
	{
	}
	
bool NxsRequiredToken::ReadValue(NxsToken &token, bool equalsAlreadyRead)
	{
	if (equalsAlreadyRead || EatEqualsThenAdvance(token))
		{
		if (EatWordThenAdvance(token, reqToken, ncl::kStringNoRespectCase))
			return WasValidRead();
		}
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a command Option with the supplied name and abbreviation
*/
NxsCmdOption::NxsCmdOption( 
  const string		  & n,			/* name of the command Option */
  bool					validityIsLabile,	/* true the test of an options validity depends on the state of the rest of the program (there is a function call to check the option) */
  bool					persist,	/* true if the Option is persistent */
  CmdPermissionLevel 	permissions) /* permissionLevel (userRestriction) that the user must reach to have access to the command */
  	:errState(no_err),
  	name(n),
  	abbreviation(n),	//the abbreviations are automatically set by the command
  	userRestriction(permissions),
  	status(at_factory), 
  	statusBefCmd(at_factory)
  	{
  	optionFlag  = (persist ? kPersistent : kNotPersitent);
  	ResetCheckIfOld(validityIsLabile);
  	}

void NxsCmdOption::ResetCheckIfOld(
  bool isLabile) /* the default value can be always valid, never valid, or might require checking */
	{
	if (isLabile)
		optionFlag |= kCheckIfOld;
	else
		optionFlag &= ~((char) kCheckIfOld);
	}
		
/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if this command Option matches the string, and the permission level sent as an argument exceeds
|	the userRestriction level placed on the command Option
*/
bool NxsCmdOption::CanReadKeyword(
  const string &s, 
  int permLevel) const
	{
	return (IsCapAbbreviation(s, GetAbbreviation()) && HasPermission(permLevel));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Must be called before a command Option is asked to parse portion of a command string. 
|	This allows the NxsCmdOption to back up it's current state (in case the command is illegal) and returns the value
|	to the default if the Option is not persistent.
|	Most derived classes will not need to override this function.
|	"decrements" the status field to the argument supplied (this status field is used to tell how old the setting that
|	a command Option holds is.
*/
void NxsCmdOption::PrepareToRead() 
	{
	statusBefCmd = status;
	if (status > old_setting)
		status = old_setting;
	ClearErrorState();
	StorePreviousValue();
	}

void NxsCmdOption::CommandWasExecuted() 
	{
	if (status > set_before_last_execute)
		{
		status = set_before_last_execute;
		if (!IsPersistent())
			ReturnToDefault();
		}
	}
	

/*----------------------------------------------------------------------------------------------------------------------
|	Typically called as the last step of reading a new value for the command Option.
|	Base-class version checks to verify that all of the Tests added to the CommandOption as requirements are met.
*/
bool NxsCmdOption::WasValidRead()
	{
	if (!parseLevelRequirements.empty())
		{
		string errM;
		if (!PerformAllTests(parseLevelRequirements, &errM))
			return FlagError(parse_req_failed, errM);
		}
	//@ need to enforce availability requirements
	return !HadError() && IsCurrentlyValid();
	}

 NxsCmdOption::~NxsCmdOption()
 	{
 	DeleteAndClearVecOfPtrs(parseLevelRequirements);
 	}

bool NxsCmdOption::ReadStringAsValue(
  string s)
  	{
  	PrepareToRead(); 
  	s << ';';
  	NxsToken token(s);
	++token;
	if (ReadValue(token, true))
		return true;
	RevertBecauseCommandFailed();
	return false;
  	}



/*----------------------------------------------------------------------------------------------------------------------
|	advances passed the = sign (if necessary) and reads the current token as the string to be stored.
|	If the command is still valid, the token stream is advance to the next token.
*/
template<>
bool SimpleCmdOption<string>::ReadValue(
 NxsToken &token,	/* the stream of tokens that are being read */
 bool equalsAlreadyRead) /* true if the equals sign has already been removed from the token stream (or is NOT expected) */
	{
	if (equalsAlreadyRead || EatEqualsThenAdvance(token))
		{
		*currentValue = token.GetTokenReference();
		//	NxsStringCmdOption will always be valid, but derived classes will want to
		//	perform validity checks, so we'll check Validity here to avoid making ReadValue virtual
		// 
		if (WasValidRead())
			{
			++token;
			return true;
			}
		}
	return false;	
	}

/*----------------------------------------------------------------------------------------------------------------------
|	advances passed the = sign (if necessary) and reads the current token as the string to be stored.
|	If the command is still valid, the token stream is advance to the next token.
|	Checks to make sure a single character is supplied as the stored value.
|	Does NOT alter the token reading to pull out single character tokens (ie this command Option expects the stored value
|	to be a single character nexus token)
*/
template<>
bool SimpleCmdOption<char>::ReadValue(
 NxsToken &token,	/* the stream of tokens that are being read */
 bool equalsAlreadyRead) /* true if the equals sign has already been removed from the token stream (or is NOT expected) */
	{
	if (equalsAlreadyRead || EatEqualsThenAdvance(token))
		{
		if (token.GetTokenLength() == 1)
			{
			*currentValue = token.GetTokenReference()[0];
			if (WasValidRead())
				{
				++token;
				return true;
				}
			}	
		else
			FlagError(unexpected_char, "single char");
		}
	return false;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Written to make output code easier to read.  Add an s to the word if plural = true OR adds  "a " before the word
|	if  includeIndefiniteArticle is true.  adds
*/
string StandardNoun(
  string retStr,
  bool includeIndefiniteArticle,
  bool plural) 
	{
	if (plural)
		{
		retStr << 's';
		return retStr;
		}
	else if (includeIndefiniteArticle)
		{
		string pref = "a ";
		pref << retStr;
		return pref;
		}
	return retStr;
	}


