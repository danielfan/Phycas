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
#include "phycas/src/ncl/nxs_defs.hpp"
#if NEW_NXS_BLOCK_AND_READER
#include "phycas/src/ncl/misc/utilities.hpp"
#include "phycas/src/ncl/nxs_block.hpp"
#include "phycas/src/ncl/nxs_token.hpp"
#include "phycas/src/ncl/nxs_exception.hpp"
#include "phycas/src/ncl/output/nxs_output.hpp"
using ncl::endl;
using std::string;

/*----------------------------------------------------------------------------------------------------------------------
| 	Checks for the semicolon that must follow the block name and then reads until END; or ENDBLOCK;
|	Uses NxsCommandManager::ProcessCommandLine() to read each command.
|	If the block implements some type of "leave" command.  It should set stopReading to true, and all subsequent commands
|	will be skipped.
|	 NOTE returning false from Read should mean that the error has been reported to the user!!
*/
bool NxsCommandManagerBlock::Read(
  NxsToken& token )	/* stream of Nexus tokens to be read */
	{
	++token;
	if(token.GetTokenReference() != ';') 
		{
		errorMsg << "Expecting ';' after " << GetID() << " block name, but found " << token.GetTokenReference() << " instead";
		throw NxsException(errorMsg, token);
		}
	for(;;)
		{
		while (token.GetTokenReference() == ';')	
			++token;
		if (token.Equals("END") || token.Equals("ENDBLOCK")) 
			{
			const string endWord(token.GetTokenReference());
			++token;
			if (token.GetTokenReference() == ';')	
				{
				NxsCommandManager::FinishedReadingBlock();
				CmdResult r = EndEncountered();
				if (r == kCmdFailedGenerateMessage)
					{
					if (errorMsg.empty())
						errorMsg << "The " << GetID() << " ended prematurely";
					throw NxsException(errorMsg, token);
					}
				return (r == kCmdSucceeded);
				}
			errorMsg = "Expecting a semicolon after ";
			errorMsg << endWord << " but found " << token.GetTokenReference() << " instead";
			throw NxsException(errorMsg, token);
			}
		else if (stopReading) //if this field is true, the rest of the commands should be ignored 
			return true;
		else 
			{
			isEmpty = false;
			try 
				{
				ProcessCommandLine(token);
				if (stopReading)
					return true;
				}
			catch (NxsX_UnknownCommand &)
				{
				if (!skipUnknownCommands)
					throw;
				SkipOneCommandLine(token);
				}
			}
		}
	}

#if defined (NCL_NXSBLOCK_SILENTLY_SKIPS_CMDS)
	/*----------------------------------------------------------------------------------------------------------------------
	|	This function is called when an unknown command named commandName is about to be skipped. This version of the 
	|	function does nothing (i.e., no warning is issued that a command was unrecognized). Override this virtual function 
	|	in a derived class to provide such warnings to the user.
	*/
	inline void NxsBlock::SkippingCommand(
	  const string &unknownCmd) 	/* the nexus token that was not recognized as a command */
		{
#		if defined(HAVE_PRAGMA_UNUSED)
#			pragma unused(unknownCmd)
#		endif
		}
#else
	/*----------------------------------------------------------------------------------------------------------------------
	|	if and outputManager has been supplied to the NxsCommandManagerBlock (in the constructor) this function will use it 
	|	called to alert the user that a command is being skipped (giving the name of the command and the block name).
	*/
	void NxsBlock::SkippingCommand(
	  const string &unknownCmd) const 	/* the nexus token that was not recognized as a command */
		{
		NxsOutputStream *outPtr = NxsOutput::GetOutputStreamPtr();
		if (outPtr != NULL)
			*outPtr << "Skipping unrecognized command (" << unknownCmd << ") in the " << GetID() << " block\n" << endl;
		}
#endif
/*----------------------------------------------------------------------------------------------------------------------
|	Called to skip to the end of a command.  Calls SkippingCommand() sending the first token from the line as an argument
|	and then reads until a semicolon (token will be left at the semicolon)
*/
void NxsCommandManagerBlock::SkipOneCommandLine(
  NxsToken &token)
  	{
  	SkippingCommand(token.GetTokenReference());
  	while (token.GetTokenReference() != ';' && !token.AtEOF())
  		++token;
  	}
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes NxsCommandManager with the OutputManager sent as an argument then calls NxsBlock::Reset
*/
NxsCommandManagerBlock::NxsCommandManagerBlock(
  const string &s,
  bool abbrev, 
  bool giveHelp, 
  bool giveSet,  
  bool giveAvailable)  
	:NxsCommandManager(abbrev, giveHelp, giveSet, giveAvailable),
	NxsBlock(s)
	{
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	sets the nexus object pointer to NULL, skipUnknownCommands to false, strictNEXUSParsing and isEnabled to true, then
|	calls NxsBlock::Reset
*/
 NxsBlock::NxsBlock(const string &s) 
  :isEnabled(true),
  skipUnknownCommands(true),
  stopReading(false),
  nexus(NULL),
  blockID(s)
	{
	ResetNxsBlockBase();
	}

#else //NEW_NXS_BLOCK_AND_READER
#include "ncl.h"

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes all pointer data members to NULL, and all bool data members to true except isUserSupplied, which is
|	initialized to false.
*/
NxsBlock::NxsBlock()
	:next(NULL),
	nexus(NULL),
	isEmpty(true),
	isEnabled(true),
	isUserSupplied(false)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version simply returns 0 but a derived class should override this function if it needs to construct
|	and run a NxsSetReader object to read a set involving characters. The NxsSetReader object may need to use this 
|	function to look up a character label encountered in the set. A class that overrides this method should return the
|	character index in the range [1..nchar].
*/
unsigned NxsBlock::CharLabelToNumber(
  string s)	/* the character label to be translated to the character's number */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(s)
#	endif
	return 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version simply returns 0, but a derived class should override this function if it needs to construct
|	and run a NxsSetReader object to read a set involving taxa. The NxsSetReader object may need to use this function to
|	look up a taxon label encountered in the set. A class that overrides this method should return the taxon index in
|	the range [1..ntax].
*/
unsigned NxsBlock::TaxonLabelToNumber(
  string s)	/* the taxon label to be translated to a taxon number */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(s)
#	endif
	return 0;
	}
#endif //NEW_NXS_BLOCK_AND_READER