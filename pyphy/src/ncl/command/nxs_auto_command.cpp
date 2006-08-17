//#include "phycas/force_include.h"
#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/command/nxs_auto_command.hpp"
#include "pyphy/src/ncl/command/nxs_cmd_param.hpp"
#include "pyphy/src/ncl/nxs_exception.hpp"
#include "pyphy/src/ncl/output/nxs_output.hpp"
using ncl::endl;
using std::string;
bool CheckCmdOptIfUnread(VecNxsCmdOptionShPtr::iterator &optIt);

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
NxsAutomaticallyParsedCommand::NxsAutomaticallyParsedCommand(
  const string &n)/* the name of the command */
  	:NxsCommand(n),
  	canReadAsterisk(false)
  	{
  	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called by ParseCommand whenever there is a fatal error during the attempt to parse a command.
|	Restores the state before the failed attempt to read a new command.
|	Shows the current settings if help was requested.
|	THROWS AN NxsException exception with the most appropriate message.
|	Checks if there is already an errormsg.  If so it is used as the NxsException message.
|	if not a default "unexpected xxx in yyy command.  Command Aborted" message is used.
*/
void NxsAutomaticallyParsedCommand::ParseFailed(
  NxsToken &tok,	/* current Token (the one that caused the error) */
  CmdResult explain)
	{
	RestoreStateBeforeRead();
	if (helpRequested)
		ShowCurrentSettings(UINT_MAX);
	if (errormsg.length() > 0)
		throw NxsException(errormsg, tok);
	errormsg << "Unexpected " << tok.GetTokenReference() << " in " << GetCurrentName() << " command.  Command Aborted.";
	throw NxsException(errormsg, tok, explain);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	if NxsCommand::skipHelpKeywords is true this function skips to the next token that is not ? or help (the command 
|	manager manipulates the NxsCommand::skipHelpKeyword)
*/
inline void	NxsAutomaticallyParsedCommand::SkipHelpRequests(
  NxsToken& token) /* the token stream */
	{
	if (NxsCommand::skipHelpKeywords)
		{
		while (token.GetTokenReference() == '?' || EqualsCaseInsensitive(token.GetTokenReference(), "HELP"))
			++token;
		}
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Called when a token is unrecognized.  If the token is too short to match a keyword unambiguously an errormsg is created
|	explaining the problem.
*/
void NxsAutomaticallyParsedCommand::CheckAbbreviationIsTooShort(
  NxsToken &token) /* the token stream */
  	{
	VecString optNames;
	for	(unsigned i= 0; i < keywords.size(); ++i)
		optNames.push_back(keywords[i]->GetName());
	VecString matches = GetVecOfPossibleAbbrevMatches(token.GetTokenReference(), optNames);
	if (matches.size() > 1)
		{
		errormsg = "Ambiguous command setting abbreviation for the ";
		errormsg << GetCurrentName() << " command:  " << token.GetTokenReference() << " matches ";
		for (unsigned i = 0; i < matches.size(); ++i)
			{
			errormsg << matches[i];
			if (matches.size() > 2 && i < matches.size() - 1)
				errormsg << ", ";
			else 
				errormsg << ' ';
			if (i == matches.size() - 2)
				errormsg << "and ";
			}
		errormsg << "keywords.\n";
		}
	}

std::string NxsCommand::GetCmdOptionNameString(
  NxsCmdOptionShPtr ci) const /* the command item that had an error */
	{
	//	compose keyword identifier string ( "XXX setting of the YYY command" )
  	//
	std::string identStr;
	bool unNamed = false;
	if (ci->GetName().length() > 0)
		{
		identStr << ci->GetName() << " setting of the ";
		unNamed = true;
		}
	identStr << GetCurrentName() << " command";
	return identStr;
	}

void NxsCommand::GetErrorFromCmdOption(
  NxsCmdOptionShPtr ci, /* the command item that had an error */
  NxsToken &tok)
	{
	errormsg = GetErrorStringFromCmdOption(ci.get(), tok, GetCmdOptionNameString(ci));
	}
/*--------------------------------------------------------------------------------------------------------------------------
|	Composes an error message by querying the command item.
*/
std::string NxsCommand::GetErrorStringFromCmdOption(
  NxsCmdOption * ci, /* the command item that had an error */
  NxsToken    & tok, /* the token stream */
  std::string   suffixStr)
  	{
  	NXS_ASSERT(ci->GetErrorState() != NxsCmdOption::no_err);
  	std::string errormsg;
	std::string transitionString;
	const int ciError = ci->GetErrorState();
	const bool unNamed = !(ci->GetName().empty());
	if (ciError == NxsCmdOption::cancelled_by_user)
		return std::string();
	// Compose error message
	// 
	if ((ciError & NxsCmdOption::miss_wd) != 0)
		{
		errormsg << "Expecting " << ci->GetErrorSnippet() << " but found " << tok.GetTokenReference();
		transitionString = " in the ";
		}
	else if ((ciError & NxsCmdOption::unrecognized) != 0 || (ciError & NxsCmdOption::unrecog_labile) != 0  || (ciError & NxsCmdOption::illegal_punctuation) != 0)
		{
		VecString validChoices = ci->GetLegalChoices();
		VecString matches = GetVecOfPossibleAbbrevMatches(tok.GetTokenReference(), validChoices);
		if (matches.size() > 1)
			{
			errormsg << "Ambiguous value for the " << ci->GetName() << " keyword:  " << tok.GetTokenReference() << " matches ";
			for (unsigned i = 0; i < matches.size(); ++i)
				{
				errormsg << matches[i];
				if (matches.size() > 2 && i < matches.size() - 1)
					errormsg << ", ";
				else 
					errormsg << ' ';
				if (i == matches.size() - 2)
					errormsg << "or ";
				}
			}
		else
			{
			if ((ciError & NxsCmdOption::illegal_punctuation) != 0)
				errormsg <<  "Illegal punctuation: the \'" << tok.GetTokenReference() << "\' character was found\n";
			else
				errormsg << "Unexpected " << tok.GetTokenReference() << "\n";
			string eSnip = ci->GetErrorSnippet();
			if (eSnip.length() == 0)
				{
				if (validChoices.size() > 0)
					errormsg << "Expecting one of the following: " << ci->GetDisplayType(true) << ' ';
				else
					errormsg << ci->GetDisplayType(true) << " was expected ";
				}
			else
				errormsg << "Expecting " << eSnip << ' ';
			}
		transitionString << (unNamed ? "in" : "after") << " the ";
		}
	else if ((ciError & NxsCmdOption::too_big) != 0 || (ciError & NxsCmdOption::too_big_labile) != 0)
		{
		errormsg << tok.GetTokenReference() << " is greater than the maximum legal value (" << ci->GetErrorSnippet() << ")";
		transitionString = "for the ";
		}
	else if ((ciError & NxsCmdOption::too_small) != 0 || (ciError & NxsCmdOption::too_small_labile) != 0)
		{
		errormsg << tok.GetTokenReference() << " is less than the minimum legal value (" << ci->GetErrorSnippet() << ")";
		transitionString= " for the ";
		}
	else if ((ciError & NxsCmdOption::read_twice) != 0)
		{
		errormsg << "The " << ci->GetName() << " setting occurred more than once";
		transitionString = " in the ";
		}
	else if ((ciError & NxsCmdOption::illegal_modulus) != 0)
		{
		errormsg << tok.GetTokenReference() << " is an illegal stride for the range " << ci->GetErrorSnippet();
		transitionString = " in the ";
		}   
	else if ((ciError & NxsCmdOption::illegal_range) != 0)
		{
		errormsg << "Illegal range, " << ci->GetErrorSnippet();
		transitionString = ", in the ";
		}		
	else if ((ciError & NxsCmdOption::reserved_word) != 0)
		{
		errormsg << tok.GetTokenReference() << " is not allowed because it is a reserved word";
		transitionString = " (error in the ";
		suffixStr << ')';
		}
	else if ((ciError & NxsCmdOption::unexpected_char) != 0)
		{
		errormsg << "Expecting a " << ci->GetErrorSnippet() << " but found " << tok.GetTokenReference();
		transitionString = " in the ";
		}
	else if ((ciError & NxsCmdOption::parse_req_failed) != 0)
		{
		errormsg  << "Could not read the value specified because " << ci->GetErrorSnippet();
		transitionString = " (error in the ";
		suffixStr << ')';
		}
	else if ((ciError & NxsCmdOption::illegal_name) != 0)
		{
		errormsg  << tok.GetTokenReference() << " is an illegal name because " << ci->GetErrorSnippet();
		transitionString = " (error in the ";
		suffixStr << ')';
		}
	else if ((ciError & NxsCmdOption::query_for_error) != 0)
		{
		errormsg << "Error:  " << ci->GetErrorSnippet();
		transitionString = " (error in the ";
		suffixStr << ')';
		}
	if (!suffixStr.empty())
		errormsg << transitionString << suffixStr;
	return errormsg << ".\n";
	}
  	
/*--------------------------------------------------------------------------------------------------------------------------
|	Tries to find a keyword for which a call of CanReadKeyword(currentToken) returns true.
|	If one is found CommandItem::TryToRead() is called and a pointer to the command item is returned.
|	if none is found, NxsCmdOptionShPtr() is returned. 
*/
NxsCmdOptionShPtr NxsAutomaticallyParsedCommand::ParseKeyword(
  NxsToken &token)
  	{
	for (VecNxsCmdOptionShPtr::iterator ci = keywords.begin(); ci != keywords.end(); ++ci)
		{
		if ((*ci)->CanReadKeyword(token.GetTokenReference(), NxsCommand::userPermLevel))
			{
			++token;
			(*ci)->TryToRead(token, false);
			return *ci;
			}
		}	
	return NxsCmdOptionShPtr();
  	}
 
bool NxsAutomaticallyParsedCommand::CheckUnnamedCmdOptIfUnread(NxsCmdOptionShPtr &cmdOpt)
	{
	const NxsCmdOption::ReadStatus &rs = cmdOpt->GetReadStatus();
	if (rs != NxsCmdOption::new_setting && OptionShouldBeValidated(cmdOpt.get()))
		{
		if ((rs == NxsCmdOption::at_factory|| cmdOpt->CheckIfOld()) && !cmdOpt->IsCurrentlyValid())
			{
			errormsg << "Error in the " << GetCurrentName() << " command.  ";
			string temp = GetUsage();
			if (temp.empty())
				{
				errormsg << "The command contains unnamed options that must are not currently documented";
				NXS_ASSERT(0);
				}
			else
				errormsg << "The correct usage of the command is:\n" << temp;				
			return false;
			}		
		}
	return true;
	}

bool NxsAutomaticallyParsedCommand::CheckNamedCmdOptIfUnread(NxsCmdOptionShPtr &cmdOpt)
	{
	const NxsCmdOption::ReadStatus &rs = cmdOpt->GetReadStatus();
	if (rs != NxsCmdOption::new_setting && OptionShouldBeValidated(cmdOpt.get()))
		{
		if (rs == NxsCmdOption::at_factory)
			{
			if (!cmdOpt->IsCurrentlyValid())
				{
				errormsg << "Error in the " << GetCurrentName() << " command.   The " << cmdOpt->GetName() << " option must be specified ";
				if (cmdOpt->IsPersistent())
					errormsg << "before the command can be executed";
				else
					errormsg << "every time the command is executed";
				return false;
				}
			}
		else if (cmdOpt->CheckIfOld() && !cmdOpt->IsCurrentlyValid())
			{
			errormsg << "Error in the " << GetCurrentName() << " command.   The setting of the " << cmdOpt->GetName() << " option is not longer valid.";	
			return false;
			}
		}
	return true;
	}

bool NxsAutomaticallyParsedCommand::AllOptionsAreValid()
	{
	VecNxsCmdOptionShPtr::iterator optIt;
	for (optIt =  unnamedCmdSettings.begin();	optIt != unnamedCmdSettings.end(); ++optIt)	
		{
		if (!CheckUnnamedCmdOptIfUnread(*optIt))
			return false;
		}		
	for (optIt =  keywords.begin();	optIt != keywords.end(); ++optIt)			
		{
		if (!CheckNamedCmdOptIfUnread(*optIt))
			return false;
		}
	return true;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Reads the token stream throwing an NxsException exception if there is an error.
|	Keywords can be in any order.
|	Unnamed settings are read in the order they were added to the command.  If one of them fails, the named keywords are used
|	to read the rest of the stream.  
|	thus if a command has keywords and 2 unnamed settings any of the following will be accepted:
|	CommandName;
|	CommandName keywords;
|	CommandName unNamed1 keywords;
|	CommandName unNamed1 unNamed2 keywords;
|	CommandName unNamed1;
|	CommandName unNamed1 unNamed2;
|	BUTNOT
|	CommandName keywords unNamed...
|	CommandName unNamed2 ... 
*/
bool	NxsAutomaticallyParsedCommand::ParseCommand(
  NxsToken& token) //throws NxsException
	{
	++token;	// token should be the command's name
	if (canReadAsterisk && token.GetTokenReference() == '*')
		{
		asteriskRead = true;
		++token;
		}
	else
		asteriskRead = false;
		
	MarkFilePosition(token);
	SkipHelpRequests(token);
	bool unNamedSuccessfullyRead = true;
	//	UnnamedSettings, if there are any, must go first, and if there are more than one,
	//	the order is specified by their order in the  unnamedCmdSettings vector (the order they
	//	were added to the command)
	//
	VecNxsCmdOptionShPtr::iterator currUnNamedSetting = unnamedCmdSettings.begin();
	for (; token.GetTokenReference() != ';' && currUnNamedSetting != unnamedCmdSettings.end(); ++currUnNamedSetting)
		{
		//	Try to read the next unnamed setting
		//
		if ((*currUnNamedSetting)->IsAvailable())
			{
			(*currUnNamedSetting)->TryToRead(token , true);	// unnamed options don't have equals signs,  so send true to indicate the equals has been read
			if ((*currUnNamedSetting)->HadError())
				{
				if ((*currUnNamedSetting)->WasErrorFatal(canExecuteHandler && ! helpRequested))
					{
					if (keywords.size() > 0)
						unNamedSuccessfullyRead = false;	// stop trying unNamed, give the named settings a shot
					else
						{
						GetErrorFromCmdOption((*currUnNamedSetting),token);
						errormsg << "Command Aborted.\n";
						ParseFailed(token, kCmdFailedGenerateMessage); //throws NxsException
						}
					}
				else
					{
					GetErrorFromCmdOption((*currUnNamedSetting),token);
					warnings.push_back(errormsg);	// non fatal error is good enough, we'll consider the token to be parsed
					errormsg.clear();
					}
				}
			if (!unNamedSuccessfullyRead)	
				break;	
			SkipHelpRequests(token);
			MarkFilePosition(token);
			}
		}
		
	if (keywords.size() > 0)
		{
		if (!unNamedSuccessfullyRead)
			{	// one of the unNamed settings had a fatal error.  
				// maybe the user did not intend any of the tokens to be unNamed.  
				// Back the token up and try to read it and reset to command item that choked
				//
			token.SeekTokenizerState(backupTokenizerState);
			NXS_ASSERT((*currUnNamedSetting)->HadError());
			(*currUnNamedSetting)->RevertBecauseCommandFailed();
			}
		SkipHelpRequests(token);
		while(token.GetTokenReference() != ';')	//walk through all of the tokens in the command
			{
			NxsCmdOptionShPtr optThatHandled = ParseKeyword(token);
			if (optThatHandled == NxsCmdOptionShPtr())
				{
				CheckAbbreviationIsTooShort(token);
				ParseFailed(token, kCmdFailedGenerateMessage); //throws NxsException
				}
			else if (optThatHandled->HadError())
				{
				GetErrorFromCmdOption(optThatHandled,token);
				if (optThatHandled->WasErrorFatal(canExecuteHandler && ! helpRequested))
					{
					errormsg << "Command Aborted.\n";
					ParseFailed(token, kCmdFailedGenerateMessage); //throws NxsException
					}
				else
					{
					warnings.push_back(errormsg);
					errormsg.clear();
					}
				}
			SkipHelpRequests(token);
			MarkFilePosition(token);
			}
		}
	if (token.GetTokenReference() != ';')
		ParseFailed(token, kCmdFailedGenerateMessage);//throws NxsException
	return true;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	writes any warning to the outputmanager (if it is not NULL)
|
*/
void NxsAutomaticallyParsedCommand::PrintWarnings() const 
	{
	NxsOutputStream * nxsWarnStream = NxsOutput::GetWarningStreamPtr();
	if (nxsWarnStream != NULL)
		{
		for (VecString::const_iterator wIt = warnings.begin(); wIt != warnings.end(); ++wIt)
			*nxsWarnStream << "Warning:\t" << *wIt << endl;
		}
	}
	
