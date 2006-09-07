//#include "phycas/force_include.h"
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/command/nxs_command_manager.hpp"
#include "phypy/src/ncl/command/nxs_cmd_param.hpp"
#include "phypy/src/ncl/command/nxs_command.hpp"
#include "phypy/src/ncl/nxs_exception.hpp"
#include "phypy/src/ncl/misc/utilities.hpp"
#include <boost/bind.hpp>
#include <boost/mem_fn.hpp>
#include "phypy/src/ncl/output/nxs_output.hpp"
#include "phypy/src/ncl/command/nxs_command_output.hpp"
#include "phypy/src/ncl/misc/string_extensions.hpp"

using std::map;
using std::vector;
using std::set;
using std::pair;
using std::string;
using ncl::endl;

void NxsCommandManager::PrintStateToHiddenStream() const
	{
	if (printStateToHiddenStream && NxsOutput::GetHiddenQueryStreamPtr())
		{
		NxsHiddenQueryStream & outStream = *NxsOutput::GetHiddenQueryStreamPtr();
		Emit<kCommandStateOutSyle>(outStream, *this) << endl;
		}
	}

void  NxsCommandManager::ProcessAliasedCommandLine(const string &expansion, NxsToken &token)
	{
	//@ the the file line, column, position will be wrong in an aliased command
	string tempStream = expansion;
	bool eofWasAllowed = token.SetEOFAllowed(false);
	while (token.GetTokenReference() != ';')
		{
		tempStream << token.GetTokenReference();
		++token;
		}
	token.SetEOFAllowed(eofWasAllowed);
	tempStream << ';';
	NxsToken fakeToken(tempStream);
	++fakeToken;
	ProcessCommandLine(fakeToken);
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
| 	Figures out which Command to pass the Command line to.
|	For normal commands: the first token is compared to the command names.  If a match is found, the NxsCommandManager calls
|	Command::CanBeRead() (to make sure hasn't been read already or can be read multiple times).  If CanBeRead returns true,
|	Command::IsAvailable() is called to see if all restrictions on the command are currently met)
|	If the command is available Command::SetAllowCommandExecution(true) and then Command::ProcessCommand() are called.
|	
|	Support for help and set commands are added here (if the commandmanager has been told to supply them)
|
|	throws NxsException exceptions if the requested command is not available.
|	XUnknownCommand (a type of NxsException exception) is thrown if the requested command is not known
|
*/
void NxsCommandManager::ProcessCommandLine(
  NxsToken	&token) /* token stream */
	{
	NXS_ASSERT(readyToParse); // If you trip this you must not have called FinishedAddingCommands
	
	if (!aliases.empty())
		{
		map<string, string, NxsStringNoCaseLess>::const_iterator aliasIt = aliases.find(token.GetTokenReference());
		if (aliasIt != aliases.end())
			{
			++token;
			ProcessAliasedCommandLine(aliasIt->second, token);
			return;
			}
		}
		
	bool helpReq = false;
	if (addHelp)
		{
		NxsCommand::skipHelpKeywords = true;
		if (PreProcessCommandLine(token, &helpReq))
			{
			if (explicitHelp)
				explicitHelp->ProcessCommand(token, false);
			else
				ProcessHelp(token);
			return;
			}
		}
	else
		NxsCommand::skipHelpKeywords = false;
		
	if (addSet && IsCapAbbreviation(token.GetTokenReference(), "SET"))
		{
		NxsTokenizerState backupState = token.GetTokenizerState();
		++token;
		//	see if the next token is a command name, if so pass the control to that command (with executeHandler set to false
		//	so that the command will be parsed but no additional action will be take
		//
		unsigned i = GetCommandIndex(token.GetTokenReference());
		if (i != UINT_MAX)
			{
			commands[i]->SetAllowCommandExecution(false);
			commands[i]->ProcessCommand(token, helpReq);
			commands[i]->SetAllowCommandExecution(true);
			return;
			}
		//	didn't see a command name.  Maybe this command manager has a different set command (with keywords other than command names)
		//	restore the state at the beginning of the command and Give the rest of this function a chance to deal with this command line
		//
		token.SeekTokenizerState(backupState);
		}
	else if (addAvailable && IsCapAbbreviation(token.GetTokenReference(), "AVAILABLE"))
		{
		ProcessAvailable(token);
		return;
		}
	//	Figure out which command was requested
	//
	unsigned i = GetCommandIndex(token.GetTokenReference());
	string eMessage;
	if (i != UINT_MAX)
		{
		//	Process all commands with the help/? in them
		//
		if (helpReq)
			{
			commands[i]->ProcessCommand(token, true);
			return;
			}
		//	make sure the command hasn't been read (if it is a once-per-block command
		//
		if (commands[i]->CanBeRead())
			{
			//	Check if the command is available (if the user level is high enough and all of the availability tests return true
			//
			if (commands[i]->IsAvailable())
				{
				commands[i]->SetAllowCommandExecution(true);
				commands[i]->ProcessCommand(token, false);
				return;
				}
			eMessage << "The " << GetCapitalized(commands[i]->GetName()) << " command is not available because " << commands[i]->GetReasonForUnavailability() << ".";
			}
		else
			eMessage << "The " << GetCapitalized(commands[i]->GetName()) << " command can only be present one time in a block.";
		throw NxsException(eMessage, token);		
		}
	//	The command is unrecognized. If we are allowing abbreviations, we should check and see if 
	//	the user got carried away and sent an abbreviation that is too short.  If so we should treat
	//	this as an error
	//	Otherwise throw and XUnknownCommand exception and the calling function determines whether
	//	unrecognized commands are errors or should just generate warnings.
	if (allowAbbreviations)
		{
		VecString commNames;
		VecString matches;
		
		for	(unsigned j= 0; j < commands.size(); ++j)
			commNames.push_back(commands[j]->GetName());
		matches = GetVecOfPossibleAbbrevMatches(token.GetTokenReference(), commNames);
		if (matches.size() > 1)
			{
			eMessage << "Ambiguous command abbreviation:  " << token.GetTokenReference() << " matches ";
			for (unsigned k = 0; k < matches.size(); ++k)
				{
				eMessage << matches[k];
				eMessage << ((matches.size() > 2 && k < matches.size() - 1) ?  ", " : " ");
				if (k == matches.size() - 2)
					eMessage << "or ";
				}
			eMessage << "commands.\n";
			throw NxsException(eMessage, token);
			}
		}
	eMessage << "The command " << token.GetTokenReference() << " is unrecognized.\n";
	throw NxsX_UnknownCommand (eMessage, token);
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	
*/
NxsCommandManager::NxsCommandManager(
  bool 	abbrev, 		/* whether or not to allow command names to be abbreviated */
  bool	giveHelp, 		/* whether or not to supply a help command */
  bool 	giveSet,		/* whether or not to supply a set command */
  bool	giveAvailable)
	:addHelp(giveHelp),
	addSet(giveSet),
	addAvailable(giveAvailable),
	allowAbbreviations(abbrev),
	readyToParse(false),
	printStateToHiddenStream(false)
	{
	}
	

/*--------------------------------------------------------------------------------------------------------------------------
| 	Adds the argument to a list of recognized commands
*/
NxsCommandShPtr NxsCommandManager::AddCommand(
  NxsCommandShPtr c)
	{
	readyToParse = false;
	NXS_ASSERT(c);
	if (c)
		{
		commands.push_back(c);
		if (StrEquals(c->GetName(), "HELP", ncl::kStringNoRespectCase))
			explicitHelp = c;
		}
	return c;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
| 	tells all of the commands that this is a new block (through Command::SetNumberOfTimesRead(0))
*/
void NxsCommandManager::PrepareToReadNewBlock()
	{
	for (unsigned i =0; i < commands.size(); ++i)
  		commands[i]->SetNumberOfTimesRead(0);
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	determines whether or not abbreviations can be used in place of the full command name.
|	returns true if all of the commands and keywords are uniquely defined
*/
bool  NxsCommandManager::SetAllowAbbreviations(
  bool letAbbrev,		/*true if the user should be allowed to abbreviate command names */
  bool passOnToCommands) /*true if you want the abbreviation status to apply to the keyword for commands (that have already been added) */
  	{
  	allowAbbreviations = letAbbrev;
  	if (passOnToCommands)
  		{
  		for (unsigned i =0; i < commands.size(); ++i)
  			commands[i]->SetAllowAbbreviatedKeywords(letAbbrev);
  		}
  	return DetermineAbbreviations();
  	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	tells all of the commands that this is a new block (through Command::SetNumberOfTimesRead(0))
*/
void NxsCommandManager::FinishedReadingBlock() const
	{
	for (unsigned i =0; i < commands.size(); ++i)
  		{
  		if (!commands[i]->WasFoundCorrectNumberOfTimes())
  			{
  			string eMessage;
  			eMessage << "The " << GetCapitalized(commands[i]->GetName()) << " command is required, but was not found.";
  			throw NxsException(eMessage);
  			}
  		}
	}
	

/*--------------------------------------------------------------------------------------------------------------------------
| 	clears list of commands (which may trigger deletion because they are stored as smart pointers)
*/
NxsCommandManager::~NxsCommandManager()
	{
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	sets the shortest unambiguous abbreviations for the commands, their items, and the choices for choice command items (if
|	they are not labile).
|	Abbreviations for some commands (e.g. help and set) may appear to be too long because they can also take any command 
|	name as settings.
*/
bool NxsCommandManager::DetermineAbbreviations()
	{
	bool retValue = true;
	if (allowAbbreviations)
		{
		vector<string> tempNames;
		//	Find the abbreviation for all the command Names
		//
		for (unsigned which_cmd = 0; which_cmd < commands.size(); ++which_cmd)
			{
			tempNames.push_back(commands[which_cmd]->GetName());
			if (!commands[which_cmd]->DetermineKeywordAbbreviations())
				retValue = false;
			}
		if (!SetToShortestAbbreviation(tempNames))
			retValue = false;
		for (unsigned i = 0; i < commands.size(); ++i)
			commands[i]->SetAbbreviation(tempNames[i]);
		}
	else
		{
		set<string> uniqueCommands;
		for (unsigned i = 0; i < commands.size(); ++i)
			{
			string capitalized = commands[i]->GetName();
			Capitalize(capitalized);
			commands[i]->SetAbbreviation(capitalized);
			pair < set<string>::iterator , bool> retVal;
			retVal = uniqueCommands.insert(capitalized);
			if (!retVal.second)
				retValue = false;
			if (!commands[i]->DetermineKeywordAbbreviations())
				retValue = false;
			}
		}
	return true;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	MUST be called after all of the command options have been added to the commands and all commands have been added to the 
|	NxsCommandManager
*/
bool NxsCommandManager::FinishedAddingCommands()
	{
	readyToParse = false;
	sort(commands.begin(), commands.end(), NamedObjLessThan<NxsCommandShPtr>());
	for (unsigned which_cmd = 0; which_cmd < commands.size(); ++which_cmd)
		{
		if (!commands[which_cmd]->FinishedAddingOptions())
			return false;
		}
	readyToParse = DetermineAbbreviations();
	NXS_ASSERT(readyToParse);	//command names weren't long enough
	return readyToParse;
	}
		
/*--------------------------------------------------------------------------------------------------------------------------
| 	returns the name of the command with index which_cmd
*/
string NxsCommandManager::GetCommandName(
  unsigned which_cmd) const
	{
	if (which_cmd >= commands.size())
		return string();	//@should throw and exception
	return commands[which_cmd]->GetName();
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	gets the index of a command by comparing the string to the commands abbreviation fields (used in ProcessCommandLine to 
|	figure out which command is being called).  UINT_MAX is returned if the string is too short to be unambiguous, or doesn't 
|	match for some other reason.
|	called for every word that could be a command in NxsCommandManager::ProcessCommandLine()
*/
unsigned	NxsCommandManager::GetCommandIndex(
  const string &n) const /* string that should match a command*/
	{
	if (allowAbbreviations)
		{
		for	(unsigned i= 0; i < commands.size(); ++i)
			{
			if (IsCapAbbreviation(n, commands[i]->GetAbbreviation()))
				return i;
			}
		}
	else 
		{
		NStrCaseInsensitiveEquals compFunc(n);
		for	(unsigned i= 0; i < commands.size(); ++i)
			{
			if (compFunc(commands[i]->GetName()))
				return i;
			}
		}
	return UINT_MAX;
	}

#if ALLOWING_MULTI_WORD_CMD_NAME
unsigned NxsCommandManager::GetCommandIndex(NxsToken &token) const
	{
	if (allowMultiWordCmdNames)
		{
		VecString prevWords;
		VecString::const_iterator pwIt;
		unsigned searchFrom;
		if (allowAbbreviations)
			{
			for	(unsigned i= 0; i < commands.size(); ++i)
				{
				searchFrom = 0;
				for (pwIt = prevWords.begin(); pwIt != prevWords.end(); ++pwIt)
#	error stopped coding here and gave up on the 	MULTI_WORD_CMD_NAMEs (temporarily)				
				if (pwIt == prevWords.end())
					{
					matchFrom token.GetTokenReference().IsCapAbbreviationFrom(
					if (IsCapAbbreviation(n, commands[i]->GetAbbreviation()))
						return i;
					}
				}
			}
		else 
			{
			NStrCaseInsensitiveEquals compFunc(n);
			for	(unsigned i= 0; i < commands.size(); ++i)
				{
				if (compFunc(commands[i]->GetName()))
					return i;
				}
			}
	
		}
	return GetCommandIndex(token.GetTokenReference());
	}
#endif
/*--------------------------------------------------------------------------------------------------------------------------
| 	Displays all of the commands in a table (if the output manager is not NULL)
*/
void NxsCommandManager::ListCommandNames()
  	{
#	if defined(NCL_HAS_TABULAR_OUTPUT)
		NxsOutputStream * outPtr = NxsOutput::GetOutputStreamPtr();
		if (outPtr != NULL)
			{
			NxsTable & outTable = *outPtr->GetTablePtr();
			outTable.Reset();
			outTable.AddString("Command");
			outTable.SetLeftMargin();
			outTable.AddString("Currently Available?");
			outTable.SetRightMargin();
			outTable.SetTopMargin();
			outTable.HyphenLine();
			for (unsigned i = 0; i < commands.size(); ++i) 
				{
				outTable.AddString(commands[i]->GetName(), 0);
				if (commands[i]->IsAvailable())
					outTable.AddString("Yes");
				else
					outTable.AddString("No");
				}
			outTable.SetBottomMargin();
			outPtr->PrintTable(&outTable);
			}
#	else
	const bool StdIOTabularOutputNotImplemented = false;
	NXS_ASSERT(StdIOTabularOutputNotImplemented); //unimplemented
#	endif
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
| 	lists command names or specific help information for the commands listed after the help request
|	Does nothing if the output manager is set to NULL.
*/
void	NxsCommandManager::ProcessHelp(
  NxsToken& token) /*command line which started with ? or help, but is now pointing to the next token */
	{
	NxsOutputStream * outPtr = NxsOutput::GetOutputStreamPtr ();
	if (outPtr != NULL)
		{
		if (token.GetTokenReference() == ';')
			ListCommandNames();
		else
			{
			bool helpForOneCmdShown = false;
			VecString unknownCmds;
			const string hStr("help");
			while (token.GetTokenReference() != ';')
				{
				unsigned i = GetCommandIndex(token.GetTokenReference());
				if (i != UINT_MAX)
					{
					helpForOneCmdShown = true;
					NxsCommandShPtr commandWithHelpReq = commands[i];
					// Check to see if the user is requesting help about a specific command keyword
					//
					++token;
					NXS_ASSERT(commandWithHelpReq);
					unsigned keywordIndex = commandWithHelpReq->GetKeywordIndex(token.GetTokenReference());
					//	display command help
					//
					if (keywordIndex != UINT_MAX)
						commandWithHelpReq->ShowDescription(keywordIndex, true);
					else
						{
						commandWithHelpReq->ShowDescription(UINT_MAX, false);
						commandWithHelpReq->ShowUsage();
						commandWithHelpReq->ShowExample(keywordIndex, (keywordIndex != UINT_MAX));
						commandWithHelpReq->ShowCurrentSettings(keywordIndex);
						}
					}
				else
					{
					if (token.GetTokenReference() == '?' || StrEquals(token.GetTokenReference(), hStr, ncl::kStringNoRespectCase))
						unknownCmds.push_back(token.GetTokenReference());
					++token;
					}
				}
			if (!unknownCmds.empty())
				{
				NxsOutputStream &nxsOutStream = *outPtr;
				nxsOutStream << "Expecting a list of command names after the help request.  The following words were unrecognized: ";
				for (VecString::iterator ucIt = unknownCmds.begin(); ucIt != unknownCmds.end(); ++ucIt)
					nxsOutStream << ' ' << *ucIt;
				ListCommandNames();
				}
			}
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	lists command names or specific help information for the commands listed after the help request
|	Does nothing if the output manager is set to NULL.
*/
void	NxsCommandManager::ProcessAvailable(
  NxsToken& token) /*command line which started with ? or help, but is now pointing to the next token */
	{
	++token;
#   if defined(NCL_PRINT_COMMAND_STATE_TO_HIDDEN) && NCL_PRINT_COMMAND_STATE_TO_HIDDEN
		PrintStateToHiddenStream();
#   else
		//deprecated available command handling
			
		NxsHiddenQueryStream * hqs = NxsOutput::GetHiddenQueryStreamPtr();
		if (hqs != NULL)
			{
			if (token.GetTokenReference() == ';')
				{
				for (unsigned i = 0; i < commands.size(); ++i) 
					{
					if (commands[i]->IsAvailable())
						*hqs << commands[i]->GetName() << '\n';
					}
				}
			else
				{
				unsigned i = GetCommandIndex(token.GetTokenReference());
				if (i != UINT_MAX)
					{
					NxsCommandShPtr commandWithHelpReq = commands[i];
					++token;
					NXS_ASSERT(commandWithHelpReq);
					VecString validArgs;
					if (token.GetTokenReference() == ';')
						{
						VecNxsCmdOptionShPtr vs = commandWithHelpReq->GetAvailableSettings();
						transform(vs.begin(), vs.end(), back_inserter(validArgs), boost::mem_fn(&NxsCmdOption::GetName));
						}
					else
						{
						unsigned keywordIndex = commandWithHelpReq->GetKeywordIndex(token.GetTokenReference());
						NxsCmdOption * cmdOpt = commandWithHelpReq->GetCmdOption(keywordIndex);
						if (cmdOpt == NULL || (++token).GetTokenReference() != ';')
							{
							string x = "Unrecognized option:  ";
							x << token.GetTokenReference();
							throw NxsException(x, token);
							}
						validArgs = cmdOpt->GetValidArgument();
						}
					for (VecString::const_iterator vaIt = validArgs.begin(); vaIt != validArgs.end(); ++vaIt) 
						*hqs << *vaIt << '\n';
					}
				}
			*hqs << ncl::endl;
		}
#   endif //1
	}

	

/*--------------------------------------------------------------------------------------------------------------------------
|	Checks for ? or "help" anywhere in the command before parsing begins, and makes sure there is a ; before the end of file
|	The token is advanced to the first non - ? and non-"help token. 
|	returns true if the first token was a help request.
*/
bool NxsCommandManager::PreProcessCommandLine(
  NxsToken &token, /* NxsToken that is reading the input stream */
  bool *helpReq)	/*reference that will be set to true if ? or help was encountered in the settings portion of the command */
	{
	string hStr = "help";
	*helpReq = false;
	
	bool isHelpComm = (token.GetTokenReference() == '?' || StrEquals(token.GetTokenReference(), hStr, ncl::kStringNoRespectCase));
	if (isHelpComm)
		{
		++token;
		if (!explicitHelp)
			{
			while (token.GetTokenReference() == '?' || StrEquals(token.GetTokenReference(), hStr, ncl::kStringNoRespectCase))
				++token;
			}
		}
	NxsTokenizerState begPos = token.GetTokenizerState();
	for (;;)
		{
		if (token.GetTokenReference() == '?' || StrEquals(token.GetTokenReference(), hStr, ncl::kStringNoRespectCase))
			*helpReq = true;
		else
			{
			if (token.GetTokenReference() == ';')
				break;
			}
		++token;
		if (token.AtEOF()) 
			{
			token.SeekTokenizerState(begPos);
			string errormsg = "Unexpected end of file encountered in ";
			errormsg << token.GetTokenReference() << " command.";
			throw NxsException(errormsg, token);
			}
		}
	token.SeekTokenizerState(begPos);
	return isHelpComm;
	}
	
