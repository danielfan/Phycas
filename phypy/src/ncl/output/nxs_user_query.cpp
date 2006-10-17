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
#include <deque>
#include <cctype> 
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/output/nxs_user_query.hpp"
#include "phypy/src/ncl/output/nxs_output_stream.hpp"
#include "phypy/src/ncl/output/nxs_sax_output_wrapper.hpp"
#include "phypy/src/ncl/output/nxs_input.hpp"	
#include "phypy/src/ncl/misc/string_extensions.hpp"
#include "phypy/src/ncl/misc/algorithm_extensions.hpp"
#include "phypy/src/ncl/command/nxs_file_cmd_param.hpp"
#include "phypy/src/ncl/command/nxs_choice_cmd_param.hpp"
#include "phypy/src/ncl/command/nxs_command.hpp"
#include "phypy/src/ncl/output/nxs_input.hpp"
#include <boost/lexical_cast.hpp>

using std::string;
using ncl::endl;
using std::list;
using NxsIO::XNxsIOError;
using NxsIO::NxsInput;
using NxsIO::UserQueryPrompter;
using boost::lexical_cast;
namespace NxsIO{
class UserQueryPrompter
	{
	public:
		virtual ~UserQueryPrompter() {}
		virtual void AlertAutoResponse(NxsWritableStream * outStream) const = 0;
		virtual void PromptUser(NxsWritableStream & outStream, const std::string & previousErr) const = 0;
		virtual bool IsCancelResponse(const std::string & ) const
			{
			return false;
			}
		virtual const std::string GetBasicPrompt() const = 0;

	};
/*			return ;
		std::string retStr = startOpenTag + " type=\"";
		switch (msgType)
			{
			case kPromptForFileChooser: 
				retStr.append("file");
				break;
			case kPromptForAnyString:
				retStr.append("string");
				break;
			case kPromptChoices:
				retStr.append("choices");
				break;
			case kPromptCancelOK:
				retStr.append("cancel_ok");
				break;
			case kPromptNoYes:
				retStr.append("no_yes");
				break;
*/	
class GeneralUserQueryPrompter: public UserQueryPrompter
	{
	public:
		GeneralUserQueryPrompter(const std::string & windowTitle, const std::string & message, const std::string & autoAcceptString)
			:_windowTitle(windowTitle),
			_message(message),
			_autoAcceptString(autoAcceptString)
			{}
		virtual ~GeneralUserQueryPrompter()
			{}
			
		virtual void AlertAutoResponse(NxsWritableStream * outStream) const
			{
			if (outStream != NULL)
				{
				*outStream << _windowTitle << '\n' << _message;
				if (!_autoAcceptString.empty())
					*outStream << '\n' << _autoAcceptString;
				*outStream << '\n' << endl;
				}
			}
		virtual void PromptUser(NxsWritableStream & outStream, const std::string & previousErr) const
			{
			const std::string promptString = previousErr + _message; 
			SetNewOutputContext(outStream, VecNxsSAXAttribute(1, NxsSAXAttribute("type", "string")));
			if (true)
				{
				NxsSaxOutputWrapper titleSOW(outStream, "title");
				titleSOW << _windowTitle;
				}
			if (true)
				{
				NxsSaxOutputWrapper msgSAXWrapper(outStream, "message");
				msgSAXWrapper << promptString ;
				}
			outStream << endl;
			}
		virtual const std::string GetBasicPrompt() const
			{
			return _message;
			}
	protected:
		std::string _windowTitle, _message, _autoAcceptString;
	};
class ChoiceUserQueryPrompter: public GeneralUserQueryPrompter
	{
	public:
		ChoiceUserQueryPrompter(const string & windowTitle, const std::string & message, const string & autoAcceptString, UInt defIndex,const NxsChoiceCmdOption & choiceOpt)
			:GeneralUserQueryPrompter(windowTitle, message, autoAcceptString),
			_defaultIndex(defIndex),
			_choiceOpt(choiceOpt)
			{
			}
		virtual ~ChoiceUserQueryPrompter()
			{}
		
		virtual void PromptUser(NxsWritableStream & outStream, const std::string & previousErr) const
			{
			const std::string promptString = previousErr + _message; 
			SetNewOutputContext(outStream, VecNxsSAXAttribute(1, NxsSAXAttribute("type", "choices")));
			if (true)
				{
				NxsSaxOutputWrapper titleSOW(outStream, "title");
				titleSOW << _windowTitle;
				}
			if (true)
				{
				NxsSaxOutputWrapper msgSAXWrapper(outStream, "message");
				msgSAXWrapper << promptString;
				}
			const VecString choices = _choiceOpt.GetLegalChoices();
			for (unsigned i = 0; i < choices.size(); i++)
				{
				VecNxsSAXAttribute atts(1, NxsSAXAttribute("index", lexical_cast<std::string>(i)));
				//VecNxsSAXAttribute atts;
				if (i == _defaultIndex)
					atts.push_back(NxsSAXAttribute("default", "true"));
				NxsSaxOutputWrapper choiceSAXWrapper(outStream, "choice", atts);
				//choiceSAXWrapper << i << ' ';
				choiceSAXWrapper << choices[i];
				}
			outStream << endl;
			}
		virtual const std::string GetBasicPrompt() const
			{
			std::string m;
			m = _message;
			const VecString choices = _choiceOpt.GetLegalChoices();
			for (unsigned i = 0; i < choices.size(); i++)
				m << '\n' << i << '\t' << choices[i];
			return m;
			}
	protected:
		const UInt   _defaultIndex;
		const NxsChoiceCmdOption & _choiceOpt;
	};
	
class FileUserQueryPrompter: public GeneralUserQueryPrompter
	{
	public:
		FileUserQueryPrompter(const std::string & message, const string & autoAcceptString)
			:GeneralUserQueryPrompter(std::string(), message, autoAcceptString)
			{}
		virtual ~FileUserQueryPrompter()
			{}
			
		virtual void AlertAutoResponse(NxsWritableStream * outStream) const
			{
			if (outStream != NULL)
				{
				*outStream << '\n' << _message;
				if (!_autoAcceptString.empty())
					*outStream << '\n' << _autoAcceptString;
				*outStream << '\n' << endl;
				}
			}
		virtual void PromptUser(NxsWritableStream & outStream, const std::string & previousErr) const
			{
			const std::string promptString = previousErr + _message; 
			SetNewOutputContext(outStream, VecNxsSAXAttribute(1, NxsSAXAttribute("type", "file")));
			if (true)
				{
				NxsSaxOutputWrapper msgSAXWrapper(outStream, "message");
				msgSAXWrapper << promptString;
				}
			outStream << endl;
			}
		virtual bool IsCancelResponse(const std::string & response) const
			{
			return EqualsCaseInsensitive(response, "CANCEL");
			}
	};
	
} //namespace NxsIO
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
/// Enters an infinite loop or prompting the user and then using a command param 
/// parser to interpret the user's response to a query. Returns if valid response
/// is read.
/// This function is used to implement the other NxsUserQuery methods.
/// \returns user's choice, or nonBlockingResponse if runWithoutInput = true
/// \throws XNxsIOError if (!runWithoutInput && !CanCommunicateWithUser()) or a NULL response is received.
/////
template<typename T, class CMD_OPT>
T   NxsUserQuery::GenericUserQuery(const UserQueryPrompter & uqp, T nonBlockingResponse, CMD_OPT & cmdOpt) const
	{
	if (runWithoutInput)
		{
		uqp.AlertAutoResponse(outStream);
		return nonBlockingResponse;
		}
	if (!CanCommunicateWithUser())
		{
		std::string s;
		s << "Error: Could not respond to the query:\n" << uqp.GetBasicPrompt() << "\nAborting.";
		throw XNxsIOError(s);
		}
	string errString;
	const std::string emptyString;
	for (;;)
		{   // append here because promptString will hold errors after the first loop through
		uqp.PromptUser(*outStream, errString);
			// get user's response
		string response;
		inputStream->AppendNextLine(&response, NULL);
		if (response.empty())
			{
			std::string s;
			s << "Error: Empty response to the query:\n" << uqp.GetBasicPrompt() << "\nAborting.";
			throw XNxsIOError(s);
			}
			// use cmdOpt to validate
		string scratch;
		cmdOpt.PrepareToRead();
		NxsToken token(response);
		++token;
			// we allow cancel responses for file types
		if (uqp.IsCancelResponse(token.GetTokenReference()))
			return nonBlockingResponse;
		cmdOpt.TryToRead(token , true);	// send true to indicate that we only want to read a value (not a <cmd param name> = <value> tuple)
		if (!cmdOpt.HadError())
			return cmdOpt.GetValue();
		errString = NxsCommand::GetErrorStringFromCmdOption(&cmdOpt, token, emptyString);
		}
	}
		

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8	
/// \returns the index of the choice that the user selected.
/// \see NxsUserQuery::GenericUserQuery
/////
unsigned NxsUserQuery::UserChoice(
  const std::string & windowTitle, /// title of alert dialog (in GUI) used to display the message
  const std::string & query, /// the query to display
  const VecString &v,  /// vector of acceptable user responses.  MUST NOT be empty
  unsigned defaultResponse,  /// the preferred response (or UINT_MAX) if there is no default MUST be less than v.size() or UINT_MAX
  unsigned nonBlockingResponse) const /// the response used in non-interactive mode to allow computation to proceed MUST be less than v.size()
	{
	
	PHYCAS_ASSERT(defaultResponse == UINT_MAX || defaultResponse < v.size());
	PHYCAS_ASSERT(nonBlockingResponse < v.size());
	string autoAcceptString;
	autoAcceptString << v[nonBlockingResponse] << " chosen  (Running in automatic mode)";
	UInt chosenIndex = defaultResponse;
	std::string scratch;
	const std::string emptyString;
	NxsChoiceCmdOption choiceCmdOpt(emptyString, &chosenIndex, &scratch, v[defaultResponse], v, true, true, kCmdPermBasicUser);
	choiceCmdOpt.SetAllowNumbers(true);
	NxsIO::ChoiceUserQueryPrompter cuqp(windowTitle, query, autoAcceptString, defaultResponse, choiceCmdOpt);
	return GenericUserQuery<unsigned int, NxsChoiceCmdOption>(cuqp, nonBlockingResponse, choiceCmdOpt);
	}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8	
/// \returns an output file path specified by the user
/// \see NxsUserQuery::GenericUserQuery
/////
NxsOutFilePath NxsUserQuery::GetOutputFilePath(
  const string & fileNameOrDescrip) const /// message to display to the user
	{
	const string autoAcceptString = "No file chosen  (Running in automatic mode)";
	NxsOutFilePath nofp;
	const std::string emptyString;
	NxsOutFileCmdOption outfCmdOp(emptyString, emptyString, &nofp, NxsOutFilePath(), NULL, true, kCmdPermBasicUser);
	string messageToDisplay;
	messageToDisplay << "Could not open \"" << fileNameOrDescrip << "\".  Select a new file or \'Cancel\'.";
	NxsIO::FileUserQueryPrompter fuqp(messageToDisplay, autoAcceptString);
	return GenericUserQuery<NxsOutFilePath, NxsOutFileCmdOption>(fuqp , NxsOutFilePath(), outfCmdOp);
	}


/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8	
/// Returns an input file path specified by the user
/////
NxsInFilePath NxsUserQuery::GetInputFilePath(
  const string & fileNameOrDescrip) const
	{
	NxsInFilePath nifp; 
	const std::string emptyString;
	NxsInFileCmdOption infCmdOp(emptyString, &nifp, NxsInFilePath(), true, kCmdPermBasicUser);
	const string autoAcceptString = "No file chosen  (Running in automatic mode)";
	string messageToDisplay;
	messageToDisplay << "Could not open \"" << fileNameOrDescrip << "\".  Select a new file or \'Cancel\'.";
	NxsIO::FileUserQueryPrompter fuqp(messageToDisplay, autoAcceptString);
	return GenericUserQuery<NxsInFilePath, NxsInFileCmdOption>(fuqp, NxsInFilePath(), infCmdOp);
	}
	
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8	
/// stores an alias to the NxsWritableStream
/////
NxsUserQuery::NxsUserQuery(
  NxsWritableStream * o, /// alias to object used to display output to the user. The NxsWritableStream MUST NOT be destroyed before the NxsUserQuery (unless SetOutputStream(NULL) is called first).
  NxsIO::NxsInput * i) /// alias to object used to read input from the user. The NxsInput MUST NOT be destroyed before the NxsUserQuery (unless SetInputStream(NULL) is called first).
  :outStream(o),
  inputStream(i),
  runWithoutInput(false)
	{
	}

// simple functions (with shorter argument lists) that indirectly call UserChoice
		

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8	
/// Displays a question to the user with "Yes" and "No" choice
/// returns true if the user clicks "Yes"
/////
bool NxsUserQuery::AskUserYesNoQuery(
  const std::string & windowTitle, /// title of alert dialog (in GUI) used to display the message
  const std::string & query) const /// the query to display
  
	{
	return 1 == UserChoice(windowTitle, query, "No|Yes", 1, 1);
	}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8	
/// Issues a warning message to the user.
/// \returns true if the user acknowledges the warning and false if the user chooses "cancel"
/// Identical to AskUserYesNo but with "OK" and "Cancel" choices
/////
bool NxsUserQuery::WarnUser(
  const std::string & windowTitle, /// title of alert dialog (in GUI) used to display the message
  const std::string & warningMessage) const /// the warning to display, should be phrased so that clicking "OK" would indicate the analysis should continue
	{
	return 1 == UserChoice(windowTitle, warningMessage, "Cancel|OK", 1, 1);
	}
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8	
/// Splits choiceString into a vector and calls other  NxsUserQuery::UserChoice().
/////
unsigned  NxsUserQuery::UserChoice(
  const std::string & windowTitle, /// title of alert dialog (in GUI) used to display the message
  const std::string & query, /// the query to display
  const std::string & choiceString, ///  a | separated list of choices to display,
  unsigned defaultResponse,  /// the preferred response (or UINT_MAX) if there is no default
  unsigned nonBlockingResponse) const /// the response used in non-interactive mode to allow computation to proceed.
	{
	VecString v = SplitString(choiceString, '|');
	return UserChoice(windowTitle, query, v, defaultResponse, nonBlockingResponse);
	}
		

