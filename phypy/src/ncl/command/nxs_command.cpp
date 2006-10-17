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
#include "phypy/src/ncl/command/nxs_command.hpp"
#include "phypy/src/ncl/misc/nxs_test.hpp" 
#include "phypy/src/ncl/command/nxs_cmd_param.hpp"
#include "phypy/src/ncl/nxs_token.hpp"
#include "phypy/src/ncl/nxs_exception.hpp"
#include "phypy/src/ncl/misc/utilities.hpp"
#include "phypy/src/ncl/output/nxs_table.hpp"
#include <boost/mem_fn.hpp>
#include "phypy/src/ncl/output/nxs_output.hpp"
using std::set;
using std::pair;
using std::string;
int NxsCommand::userPermLevel = kCmdPermBasicUser;		
bool NxsCommand::skipHelpKeywords = false;

VecConstNxsCmdOptionPtr NxsCommand::GetAllSettings() const
	{
	VecNxsCmdOptionShPtr::const_iterator sIt;
	VecConstNxsCmdOptionPtr v;
	v.reserve(unnamedCmdSettings.size() + keywords.size());
	for (sIt = unnamedCmdSettings.begin(); sIt != unnamedCmdSettings.end(); ++sIt)
		v.push_back(sIt->get());
	for (sIt = keywords.begin(); sIt != keywords.end(); ++sIt)
		v.push_back(sIt->get());
	return v;
	}

NxsCmdOptionShPtr NxsCommand::GetEqualsSignRequiredTokenShPtr()
	{
	return NxsCmdOptionShPtr(new NxsRequiredToken("="));
	}

/*---------------------------------------------------------------------------------------------------
|	Calls NxsCommand::PrepareToRead() to back up the state of command (in case there is an error)
|	calls ParseCommand().  If no exceptions is thrown. the setting requirements are checked.  If they
|	are passed ( and the NxsCommandManager hasn't told the command not to execute) then the execution
|	requirements are checked.  If the execution requirements are passed the ExecuteControlledFunction()
|	is called.  If it returns true, the command notes that the last read was successful.
|	Any failure in this process will result in an NxsException exception being thrown.
*/
void NxsCommand::ProcessCommand(
  NxsToken &token,
  bool helpReq) /* whether or not the user is requesting help (if true,  ? and help tokens should be skipped and instead of executing the controlled function, help should be shown) */
	{
	PrepareToRead();
	lastReadWasSuccessfulExecute = false;
	helpRequested = helpReq;
	try
		{
		ParseCommand(token); // pure-virtual
		std::string errM;
		if (!PerformAllTests(setRequirements, &errM))
			{
			errormsg = ComposeErrorPrefix();
			errormsg << errM;
			AbortCommand();	
			}
		if (helpRequested)
			{
			PrintWarnings();
			ShowCurrentSettings(UINT_MAX);
			}	
		else if (canExecuteHandler)
			{
			if (!PerformAllTests(executeRequirements, &errM))
				{
				errormsg = ComposeErrorPrefix();
				errormsg << errM;
				AbortCommand();	
				}
			if (!AllOptionsAreValid())
				AbortCommand();	
			CmdResult cmdResCode = ExecuteControlledFunction();
			if (cmdResCode == kCmdSucceeded) // pure-virtual 
				{
				++numTimesRead;
				lastReadWasSuccessfulExecute = true;
				}
			else
				AbortCommand(cmdResCode);
			}
		else
			PrintWarnings();
		}
	catch (...)
		{
		RestoreStateBeforeRead();
		throw;
		}
	}

/*---------------------------------------------------------------------------------------------------
|	Creates a command with the supplied name which is available to the basic_user level.
*/
NxsCommand::NxsCommand(
  const string &n)
  	:origCmdName(n),
  	asteriskRead(false),
  	allowAbbrevKeywords(false),
  	canExecuteHandler(true),
  	helpRequested(false),
  	lastReadWasSuccessfulExecute(true),
  	backupTokenizerState(),
  	numTimesRead(0),
  	numTimesCmdExpected(kNoLimit),
  	requiredUserLevel(kCmdPermBasicUser)
  	{
	}

NxsCommand::~NxsCommand()
	{
	DeleteAndClearVecOfPtrs(availableRequirements);
	DeleteAndClearVecOfPtrs(setRequirements);
	DeleteAndClearVecOfPtrs(executeRequirements);
	DeleteAndClearVecOfPtrs(keywords);
	DeleteAndClearVecOfPtrs(unnamedCmdSettings);
	}
	
/*---------------------------------------------------------------------------------------------------
|	Writes the header for the current settings table shown in response to help requests
*/
void NxsCommand::DisplayCurrentSettingHeader(
  NxsTable &outTable)
	{
	outTable.AddString("Setting"); // let's not invent any command names longer than max_command_len!
	outTable.SetLeftMargin();
	outTable.AddString("Current");
	outTable.AddString("Abbr.", 0); // 0 for width means make column as narrow as possible
	outTable.AddString("Type"); // setting type_colw as maximum column width, even if string is longer
	outTable.SetRightMargin();
	outTable.HyphenLine();
	outTable.SetTopMargin();
	}
/*---------------------------------------------------------------------------------------------------
|	Writes the current settings into the table shown in response to help requests
*/
void NxsCommand::DisplayKeywordCurrentSetting(
  NxsTable &outTable, 
  NxsCmdOptionShPtr ci)
	{
	outTable.AddString(ci->GetName());
	outTable.AddString(ci->GetCurrentValueAsString());
	string temp = UpperCasePrefix(ci->GetAbbreviation());
	ToLower(temp);
	outTable.AddString(temp, 0);
	outTable.AddString(ci->GetDisplayType());
	}
	
/*---------------------------------------------------------------------------------------------------
|	Backs up the command items and the command to the previous state
*/
void NxsCommand::RestoreStateBeforeRead()
	{
	VecNxsCmdOptionShPtr::iterator ci;
	for (ci = unnamedCmdSettings.begin(); ci != unnamedCmdSettings.end(); ++ci)
		(*ci)->RevertBecauseCommandFailed();
	for (ci = keywords.begin(); ci != keywords.end(); ++ci)
		(*ci)->RevertBecauseCommandFailed();
	helpRequested = false;
	}
	
/*---------------------------------------------------------------------------------------------------
|	Informs the command item that a command is about to be read (so they can return to their default 
|	settings if they aren't persistent.
*/
void NxsCommand::PrepareToRead()	
	{
	errormsg.clear();
	for_each(unnamedCmdSettings.begin(), unnamedCmdSettings.end(), boost::mem_fn(&NxsCmdOption::PrepareToRead));
	for_each(keywords.begin(), keywords.end(), boost::mem_fn(&NxsCmdOption::PrepareToRead));
	}

void NxsCommand::NotifyCmdOptionsOfExecution()
	{
	// first we check if the handler modified the command (currently this is only true for situations in which the user
	//	is prompted for replace/append for a file.  The OutFileCmdOption detects the change and returns the correct
	//	settings from GetHandlersChangesToCmd()
	//
	string r;
	for (VecNxsCmdOptionShPtr::iterator kIt = keywords.begin(); kIt != keywords.end(); ++kIt)
		r << (*kIt)->GetHandlersChangesToCmd();
	if (!r.empty())
		{
		// Read the new changes in "set" mode
		r << ';';
		NxsToken token(r);
		SetAllowCommandExecution(false);
		ProcessCommand(token, false);
		SetAllowCommandExecution(true);
		}
	
	//  this is where we need to add a hook for gathering a string that summarizes the command's state at the execution 
	// (for logging commands issued)
	for_each(unnamedCmdSettings.begin(), unnamedCmdSettings.end(), boost::mem_fn(&NxsCmdOption::CommandWasExecuted));
	for_each(keywords.begin(), keywords.end(), boost::mem_fn(&NxsCmdOption::CommandWasExecuted));
	}
	
/*---------------------------------------------------------------------------------------------------
|	returns the error string "CurrentCommandName aborted because "
*/
string NxsCommand::ComposeErrorPrefix()
	{
	string s = GetCurrentName();
	s << " aborted because ";
	return s;
	}
		
/*--------------------------------------------------------------------------------------------------------------------------
|	Called by ProcessCommand whenever there is a fatal error in the attempt AFTER the command was parsed.  
|	Restores the state before the failed attempt to read a new command.
|	Shows the current settings if help was requested.
|	THROWS AN NxsException exception with the most appropriate message.
|	ASSUMES that there is already an errormsg.  
*/
void NxsCommand::AbortCommand(
  CmdResult cResCode)	/*whether or not to explain the failure of the command (false implies that the user has already been alerted to the problem)*/
	{
	if (cResCode == kCmdFailedGenerateMessage)
		{
		PrintWarnings();
		NXS_ASSERT (errormsg.length() > 0);
		}
	else
		errormsg.clear();
	RestoreStateBeforeRead();
	NxsTokenPosInfo posInfo = backupTokenizerState.GetPosInfo();
	throw NxsException(errormsg, posInfo.GetFilePosition(), posInfo.GetFileLine(), posInfo.GetFileColumn(), cResCode);
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Uses the NxsOutputStream  (if it is not NULL) to show all of the known examples of the commands usage.
*/
void NxsCommand::ShowExample(
  unsigned , /* not supported should be the command option that is included in the example */ 
  bool requestedExplicitly) const
	{
	NxsOutputStream * outStreamPtr = NxsOutput::GetOutputStreamPtr();
	if (outStreamPtr != NULL)
		{
		NxsOutputStream &nxsOutStream = *outStreamPtr;
		if (examples.size() > 0)
			{
			nxsOutStream << "\nExample:\n\n";
			for (unsigned i = 0 ; i < examples.size(); ++i)
				nxsOutStream << examples[i];
			}
		else if (requestedExplicitly)
			nxsOutStream << "\n\n" << "Sorry, there are no documented examples of the " << GetName() << " command.";
		nxsOutStream << ncl::endl;
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Uses the NxsOutputStream (if it is not NULL) to show the description of the command
*/
void NxsCommand::ShowDescription(
  unsigned keywordIndex, 
  bool requestedExplicitly) const
	{
	NxsOutputStream * outStreamPtr = NxsOutput::GetOutputStreamPtr();
	if (outStreamPtr != NULL)
		{
		NxsOutputStream &nxsOutStream = *outStreamPtr;
		if (keywordIndex != UINT_MAX)
			{
			NXS_ASSERT(keywordIndex < (unsigned) keywords.size());
			NxsCmdOptionShPtr argCI = keywords[keywordIndex];
			if (argCI->GetDescription().length() > 0)
				nxsOutStream << argCI->GetDescription() << '\n';
			else if (requestedExplicitly)
				nxsOutStream << "\n\nSorry, the " << argCI->GetName() << " of the " << GetName() << " command is not documented." << "\n";
			}
		else 
			{
			if (GetDescription().length() > 0)
				nxsOutStream << "\n\n" << GetDescription() << "\n";
			else if (requestedExplicitly)
				nxsOutStream << "\n\nSorry, the " << GetName() << " command is not described." << "\n";
			}
		nxsOutStream << ncl::endl;
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Uses the NxsOutputStream (if it is not NULL) to show the current settings of the command in a table
*/
void NxsCommand::ShowCurrentSettings(
  unsigned keywordIndex)
  	{
	NxsOutputStream * outStreamPtr = NxsOutput::GetOutputStreamPtr();
	if (outStreamPtr != NULL)
		{
#		if defined(NCL_HAS_TABULAR_OUTPUT)
			NxsOutputStream &nxsOutStream = *outStreamPtr;
			NxsTable &outTable = *nxsOutStream.GetTablePtr();
			if (keywords.size() > 0)
				{
				// return non persistent commands to their defaults
				PrepareToRead();
				nxsOutStream << "\nCurrent settings of the " << GetName();
				if (UpperCasePrefix(abbrevCmdName).length() < abbrevCmdName.length())
					nxsOutStream << " (" << UpperCasePrefix(abbrevCmdName) << ") ";
				nxsOutStream << "command:\n";
					
				outTable.Reset();
				DisplayCurrentSettingHeader(outTable);
				if (keywordIndex == UINT_MAX)
					{
					VecNxsCmdOptionShPtr::iterator ci;
					for (ci = keywords.begin(); ci != keywords.end(); ++ci)
						{
						if ((*ci)->HasPermission(userPermLevel))
							DisplayKeywordCurrentSetting(outTable,*ci);
						}
					}
				else
					{
					NXS_ASSERT (keywordIndex < keywords.size());
					DisplayKeywordCurrentSetting(outTable, keywords[keywordIndex]);
					if (!keywords[keywordIndex]->HasPermission(userPermLevel))
						{
						nxsOutStream << "That option is only available when the user level is set to ";
						if (keywords[keywordIndex]->GetPermissionLevel() == kCmdPermAdvancedUser)
							nxsOutStream << "advanced.\n";
						else if (keywords[keywordIndex]->GetPermissionLevel() == kCmdPermDeveloper)
							nxsOutStream << "developer.\n";
						}
					}
				outTable.SetBottomMargin();
				nxsOutStream.PrintTable(&outTable);
				}
			else
				{
				nxsOutStream << "\nThere are no optional settings for the " << GetName();
				if (UpperCasePrefix(abbrevCmdName).length() < abbrevCmdName.length())
					nxsOutStream << " (" << UpperCasePrefix(abbrevCmdName) << ") ";
				nxsOutStream << "command.\n";
				}
			if (!IsAvailable())
				nxsOutStream << "This command is currently unavailable because " << GetReasonForUnavailability() << '.';
			nxsOutStream << ncl::endl;
#		else
			NXS_ASSERT(0); //@unimplemented
#		endif		
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Uses the NxsOutputStream (if it is not NULL) to show the current settings of the command in a table
*/
VecNxsCmdOptionShPtr NxsCommand::GetAvailableSettings() 
	{
	PrepareToRead();
	VecNxsCmdOptionShPtr retVec;
	copy_if(keywords.begin(), keywords.end(), std::back_inserter(retVec),  boost::mem_fn(&NxsCmdOption::IsAvailable));
	return retVec;
	}

void NxsCommand::ShowUsage() const
	{
	const string temp = GetUsage();
	if (!temp.empty())
		NxsOutput::GetNullCheckingStream<NxsOutput::kOutput>() <<  "\nUsage:\n" << temp << ncl::endl;
	}

string NxsCommand::GetUsage() const
	{
	string retStr;
	if (usage.empty())
		{
		if (unnamedCmdSettings.size() > 0)
			{
			//	There are unnamedCommandOptions in this command, if they all have descriptions, display the descriptions
			//
			string msg;
			for (unsigned i = 0; i < unnamedCmdSettings.size(); ++i)
				{
				if (unnamedCmdSettings[i]->GetDescription().length() > 0)
					msg << " <" << unnamedCmdSettings[i]->GetDescription() << '>';
				else	
					return retStr;
				}
			retStr << GetName()  << msg << ';';
			}
		// if there are no unnamed options the usage should be obvious from help command-name 
		//	(if it isn't obvious, a usage string should have been supplied when creating this command)
		//
		return retStr;
		}
	else
		retStr << GetName() << ' ' << usage;
	return retStr;
	}	
/*--------------------------------------------------------------------------------------------------------------------------
|	returns the index of the command s, or -1 if the string doesn't match any of the commands
*/
unsigned NxsCommand::GetKeywordIndex(
  const string &s) const
	{
	if (allowAbbrevKeywords)
		{
		for	(unsigned i= 0; i < keywords.size(); ++i)
			{
			if (IsCapAbbreviation(s, keywords[i]->GetAbbreviation()))
				return i;
			}
		}
	else
		{
		for	(unsigned i= 0; i < keywords.size(); ++i)
			{
			if (EqualsCaseInsensitive(s, keywords[i]->GetAbbreviation()))
				return i;
			}
		}
	
	return UINT_MAX;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	called after all options have been added to determine minimal abbreviations
*/
bool NxsCommand::FinishedAddingOptions()
	{
	sort(keywords.begin(), keywords.end(), NamedObjLessThan<NxsCmdOptionShPtr>());
	return DetermineKeywordAbbreviations();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Determines the shortest unambiguous abbreviations of keywords and choices in multiple-choice style settings
*/
bool NxsCommand::DetermineKeywordAbbreviations()
	{
	bool retValue = true;
	VecNxsCmdOptionShPtr::iterator ci;
	if (allowAbbrevKeywords)
		{
		VecString tempNames;
		for (ci = keywords.begin(); ci != keywords.end(); ++ci)
			{
			tempNames.push_back( (*ci)->GetName());
			string convN = (*ci)->GetConverseName();
			if (convN.length() > 0)
				tempNames.push_back( convN );
			if (!(*ci)->AbbreviateChoices())
				retValue = false;
			}
		if (!SetToShortestAbbreviation(tempNames))
			retValue = false;
		// now put the abbreviations back in the same order
		//
		VecString::iterator tnIt = tempNames.begin();
		for (ci = keywords.begin(); ci != keywords.end(); ++ci)
			{
			(*ci)->SetAbbreviation(*tnIt);
			++tnIt;
			// remember that bools put in 2 names so they need to take out 2 names
			//
			if ((*ci)->GetConverseName().length() > 0)
				{
				(*ci)->SetConverseAbbreviation(*tnIt);
				++tnIt;
				}
			}
		NXS_ASSERT(tnIt == tempNames.end());
		}
	else
		{
		set<string> uniqueKeywords;
		for (ci = keywords.begin(); ci != keywords.end(); ++ci)
			{
			string tempStr = (*ci)->GetName();
			Capitalize(tempStr);
			pair < set<string>::iterator , bool> retVal;
			retVal = uniqueKeywords.insert(tempStr);
			if (!retVal.second)
				retValue = false;
			(*ci)->SetAbbreviation(tempStr);
			tempStr = (*ci)->GetConverseName();
			if (tempStr.length() > 0)
				{
				Capitalize(tempStr);
				retVal = uniqueKeywords.insert(tempStr);
				if (!retVal.second)
					retValue = false;
				(*ci)->SetConverseAbbreviation(tempStr);
				}
			}
		}
	return retValue;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Displays the usage string or generates a default command usage message (and saves it as the usage string)
*/
/*--------------------------------------------------------------------------------------------------------------------------
|	returns "it is restricted to {advanced users|registered users|beta testers|developers}" depending on the argument.
*/
string UserLevelTooLowString(int requiredUserLevel)
 	{
 	string r("it is restricted to ");
	if (requiredUserLevel == kCmdPermAdvancedUser)
		r << "advanced users";
	else if (requiredUserLevel == kCmdPermRegisteredUsers)
		r << "registered users";
	else if (requiredUserLevel == kCmdPermBetaTesters)
		r << "beta testers";
	else
		r << "developers";
	return r;
	}
		
NxsCmdOptionShPtr  NxsCommand::GetKeyword(const string &n)
	{
	for (VecNxsCmdOptionShPtr::const_iterator kIt = keywords.begin(); kIt != keywords.end(); ++kIt)
		{
		if (EqualsCaseInsensitive(n, (*kIt)->GetName()))
			return *kIt;
		}
	return NxsCmdOptionShPtr();	
	}

/*---------------------------------------------------------------------------------------------------
|	returns the command name (preceded by set or help if applicable)
*/
string NxsCommand::GetCurrentName() const
	{
	string s;
	if (!canExecuteHandler)
		s << (helpRequested ? "Help for ": "Set ");
	s << GetName();
	return s;
	}

/*---------------------------------------------------------------------------------------------------
|	returns true if the user's permission level is high enough and all availablity requirements for 
|	the command are passed.
*/
bool NxsCommand::IsAvailable() const
	{
	if (userPermLevel < requiredUserLevel)
		return false;
	return PerformAllTests(availableRequirements, NULL);
	}

/*---------------------------------------------------------------------------------------------------
|	returns true if the user's permission level is high enough and all availablity requirements for 
|	the command are passed.
*/
std::string NxsCommand::GetReasonForUnavailability() const
	{
	if (userPermLevel < requiredUserLevel)
		return UserLevelTooLowString(requiredUserLevel);
	std::string errM;
	PerformAllTests(availableRequirements, & errM);
	return errM;
	}


