#ifndef NCL_NXS_COMMAND_H
#define NCL_NXS_COMMAND_H

#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/misc/algorithm_extensions.hpp"
#include "phypy/src/ncl/command/nxs_command_decl.hpp"
#include "phypy/src/ncl/nxs_token.hpp"
class NxsTable;

/*---------------------------------------------------------------------------------------------------
|	Abstract base class.  Commands should be one of the three derived classes ManualCommand, 
|	AutoCommand, or AutoSettingCommand
|	Parses commands and calls the appropriate function if the command is legal.
|	The NxsCommandManager calls Command::ProcessCommand() after it encounters the commands name.
|	ParseCommand and ExecuteControlledFunction pure virtual functions which are overridden in 
|	ManualCommand and AutoCommand.  
|	
|	The general model for a command processing is ParseCommand(NxsToken) which is expected to read
| 	the entire command line.  If there are no errors a "controlled function" is then called.
|	The controlled function performs the action that the user requested, but doesn't need to be 
|	cluttered by code checking that the values for settings are legal (ParseCommand should guarantee
|	this)
|
|	// to do : design the generic interface between commands and their GUI counterparts (menu items
|	or dialog boxes)
|
|	Note Tests and CmdOptions are sent as boost::shared_ptrs
*/
class NxsCommand
	{
	public :
		enum 	NumTimesPerBlock 	// arguments to SetNumberOfTimesAllowedPerBlock() command
				{ 	
				kOnceAlways, 	// every block must have this command 
				kOnceMax, 		// the command can occur 0 or 1 times
				kAtLeastOnce,	// command must be present (may be encountered several times
				kNoLimit			
				};	
				
		NxsCommand(const std::string &);
		virtual ~NxsCommand();
		
		// Initialization
		//
		void					AddCommandAvailabilityReq(NxsTestShPtr cr);
		void					AddSettingReq(NxsTestShPtr cr);
		void					AddExecutionReq(NxsTestShPtr cr);
		void					AddKeyword(NxsCmdOptionShPtr k, bool alwaysVerify = true);
		void					AddUnnamedSetting(NxsCmdOptionShPtr u, bool alwaysVerify = true);
		bool 					AsteriskRead()	const {return asteriskRead;}
		void					ExpectEqualsSign();
		bool					FinishedAddingOptions();
		void					SetAllowAbbreviatedKeywords(bool letAbbreviate);
		inline void				SetDescription(const std::string &);
		
		//	Accessors
		//
		bool					CanBeRead() const;
		const std::string	  & GetAbbreviation() const;
		VecConstNxsCmdOptionPtr GetAllSettings() const; 
		VecNxsCmdOptionShPtr	GetAvailableSettings(); 
		const NxsCmdOption    * GetCmdOption(unsigned keywordIndex) const
			{
			return _GetCmdOptionImpl(keywordIndex);
			}
		
		NxsCmdOption * GetCmdOption(unsigned keywordIndex)
			{
			return _GetCmdOptionImpl(keywordIndex);
			}
			
		std::string					GetCurrentName() const;
		inline const std::string  & GetDescription() const;
		std::string					GetErrorMessage() const 
									{
									return errormsg;
									}
		NxsCmdOptionShPtr		GetKeyword(const std::string &n);
		const std::string     & GetName() const;
		NumTimesPerBlock GetNumberOfTimesAllowedPerBlock() const	
			{
			return numTimesCmdExpected;
			}
		int GetNumberOfTimesRead() const	
			{
			return numTimesRead;
			}
		std::string				GetReasonForUnavailability() const;
		const VecNxsCmdOptionShPtr & GetKeywords() const  
			{
			return keywords;
			}
		const VecNxsCmdOptionShPtr & GetUnnamedSettings() const  
			{
			return unnamedCmdSettings;
			}
		std::string				GetUsage() const;
		bool					OptionShouldBeValidated(const NxsCmdOption *) const;
		bool					IsAvailable() const;
		bool					WasFoundCorrectNumberOfTimes() const;

		// Error reporting
		STATELESS_FUNC std::string GetErrorStringFromCmdOption(NxsCmdOption *, NxsToken &, std::string identString);
		void					GetErrorFromCmdOption(NxsCmdOptionShPtr, NxsToken &);
		std::string				GetCmdOptionNameString(NxsCmdOptionShPtr) const;

		// Modifiers
		//
		void					SetAllowCommandExecution(bool canExecute = true);
		void					SetNumberOfTimesAllowedPerBlock(NumTimesPerBlock neor)	{numTimesCmdExpected = neor;}
		void					SetNumberOfTimesRead(int ntr)	{numTimesRead = ntr;}
		void 					SetUsage(const std::string &);
		
		virtual void 			ShowCurrentSettings(unsigned keywordIndex);
		virtual void			ShowDescription(unsigned keywordIndex, bool reqExplicitly) const;
		virtual void			ShowExample(unsigned keywordIndex, bool reqExplicitly) const;
		virtual void			ShowUsage() const;
		//	Utilities
		//
		void 					ProcessCommand(NxsToken &token, bool helpReq);
		
		STATIC_DATA int			userPermLevel;		/* current permissions level (value from the perms enum in NxsCmdOption).  Determines which options are available */
	protected :
		STATIC_DATA bool		skipHelpKeywords;

		//	Pure virtuals (see AutoCommand and ManualCommand)
		//
		virtual CmdResult		ExecuteControlledFunction() = 0;
		virtual bool			ParseCommand(NxsToken&) = 0;

		virtual bool			AllOptionsAreValid() {return true;} //called right before ExecuteControlledFunction.  provides a hook to make sure that all of the required settings were listed in the command
		void					AbortCommand(CmdResult explain = kCmdFailedGenerateMessage);	//throws NxsException
		
		std::string				ComposeErrorPrefix();
		bool					DetermineKeywordAbbreviations();
		unsigned				GetKeywordIndex(const std::string &s) const;	
		
		void					MarkFilePosition(const NxsToken &currToken);
		void					NotifyCmdOptionsOfExecution();
	
		virtual void 			PrepareToRead();
		virtual	void 			PrintWarnings() const; // base class version does nothing
		void					RestoreStateBeforeRead();
		void					SetAbbreviation(const std::string &n);
		
		//	Displaying Help
		//	
		void					DisplayCurrentSettingHeader(NxsTable &);
		void					DisplayKeywordCurrentSetting(NxsTable &, NxsCmdOptionShPtr );
		
		STATELESS_FUNC NxsCmdOptionShPtr	GetEqualsSignRequiredTokenShPtr();
		
		typedef std::vector<NxsTestShPtr> VecTestShPtr;
		
		std::string						origCmdName;			/* the command name with capitalization as specified in HandlerObject::SetupCommands() (used for display to the user) */
		std::string						abbrevCmdName;			/* the command name with upper case prefix indicating the minimum necessary portion */
		VecNxsCmdOptionShPtr			keywords;				/* list of keywords recognized by the command */
		VecNxsCmdOptionShPtr			unnamedCmdSettings;		/* list of unnamed settings that should follow the command's name */
		std::set<const NxsCmdOption *>	unverifiedSettings;
		
		bool					asteriskRead;			/* true if an asterisk was read as the first token (only applicable for some commands see AutoCommand::SetReadsAsterisk) */
		bool					allowAbbrevKeywords;	/* whether or not abbreviations of command keywords are accepted*/
		bool					canExecuteHandler;		/* true if the command being read is not set or help.  If the command has no errors and executeHandler is true, the handler function will be called */
		bool					helpRequested;			/* true if there was a ? in the command line.  Indicates : don't execute the handler AND be a little more lenient in parsing the command (don't throw abort without displaying help)*/
		bool					lastReadWasSuccessfulExecute; /* true if the previous reading of the command resulted in exectution of the controlled function (if true, non persistent options should be reset)*/
		NxsTokenizerState		backupTokenizerState;	/* location in the token stream when MarkFilePosition was last called*/
		int						numTimesRead;			/* stores the number of times the command has been successfully read and executed since SetNumberOfTimesRead(0) was called */
		//	requirements for completion of the command
		//
		NumTimesPerBlock			numTimesCmdExpected;   	/* the number of times the command should appear in each block */
		int							requiredUserLevel;		/* the level from the NxsCmdOption permissions enum that the user must be at to have access to this command */
		mutable VecTestShPtr		availableRequirements;	/* a collection of tests that must be passed before the handlerFunction is called and do NOT depend on the settings of the command (only test whether the command is available given the current state of the program)*/
		VecTestShPtr 				setRequirements;		/* a collection of tests that must be passed before the handlerFunction is called and do NOT depend on the settings of the command (only test whether the command is available given the current state of the program)*/
		VecTestShPtr 				executeRequirements;	/* a collection of tests that must be passed before the handlerFunction is called and do NOT depend on the settings of the command (only test whether the command is available given the current state of the program)*/

		//	error reporting/help
		//
		std::string 					errormsg;				/* string describing errors associated with the current attempt to read the command */
		std::string 					descrip;				/* string describing to the user what the command does */
		std::string 					usage;					/* for odd commands this field stores a string explaining the usage */
		VecString					examples;				/* a vector of string giving examples of the command's usage */
		
		friend class NxsCommandManager;
	private :
		void				AddParamImpl(VecNxsCmdOptionShPtr * v, NxsCmdOptionShPtr p, bool alwaysVerify);
		NxsCommand(const NxsCommand &);			//never use	- don't define - aliasing makes simple element copy incorrect
		NxsCommand &operator=(const NxsCommand&);	//never use	- don't define - aliasing makes simple element copy incorrect
		NxsCmdOption * _GetCmdOptionImpl(unsigned keywordIndex) const
			{
			if (keywordIndex >= keywords.size())
				return NULL;
			return keywords[keywordIndex].get();
			}	
	};
	

std::string UserLevelTooLowString(int);

inline void NxsCommand::SetUsage(
  const std::string &u)
	{
	usage = u;
	}
		
/*---------------------------------------------------------------------------------------------------
|	Sets the backupTokenizerState to the current position in the token stream.
*/
inline void NxsCommand::MarkFilePosition(
  const NxsToken &token)
  	{
  	backupTokenizerState = token.GetTokenizerState();
  	}
  
/*---------------------------------------------------------------------------------------------------
|	Adds a test that must be passed for the command to be available for reading
*/
inline void NxsCommand::AddCommandAvailabilityReq(NxsTestShPtr cr)  
	{
	availableRequirements.push_back(cr);
	}
	
/*---------------------------------------------------------------------------------------------------
|	Adds a test that must be passed for the command to be have any setting changed
*/
inline void NxsCommand::AddSettingReq(NxsTestShPtr cr) 
	{
	setRequirements.push_back(cr);
	}
	
/*---------------------------------------------------------------------------------------------------
|	Adds a test that must be passed for the command's controlled function to be executed.
*/
inline void NxsCommand::AddExecutionReq(NxsTestShPtr cr) 
	{
	executeRequirements.push_back(cr);
	}
	
/*---------------------------------------------------------------------------------------------------
|	Adds a named command item (keyword) that can follow the command.
|	Keywords can occur in any order in the user's command string.  They are expected to be
|	in the form : keyword = settingValue
*/
inline void NxsCommand::AddKeyword(NxsCmdOptionShPtr k, bool alwaysVerify)	
	{
	AddParamImpl(&keywords, k, alwaysVerify);
	}

inline bool NxsCommand::OptionShouldBeValidated(const NxsCmdOption * c) const
	{
	return (unverifiedSettings.find(c) == unverifiedSettings.end());
	}
/*---------------------------------------------------------------------------------------------------
|	Adds a unnamed command item.  The order that unnamed command items are added to the command 
|	determines the order they must appear.
|	an example is CharSet <IndexSetItem>;
|	Note no = sign is expected before the command item
*/
inline void NxsCommand::AddUnnamedSetting(NxsCmdOptionShPtr u, bool alwaysVerify)	
	{
	AddParamImpl(&unnamedCmdSettings, u, alwaysVerify);
	}

/*---------------------------------------------------------------------------------------------------
|	Adds param
*/
inline void NxsCommand::AddParamImpl(VecNxsCmdOptionShPtr * v, NxsCmdOptionShPtr p, bool alwaysVerify)	
	{
	v->push_back(p);
	if (!alwaysVerify)
		unverifiedSettings.insert(p.get());
	}
	
/*---------------------------------------------------------------------------------------------------
|	Returns the command name with the minimal abbreviation capitalized
*/
inline const std::string & NxsCommand::GetAbbreviation() const 
	{
	return abbrevCmdName;
	}
	
/*---------------------------------------------------------------------------------------------------
|	Returns the command name (as specified in the constructor)
*/
inline const std::string & NxsCommand::GetName() const	
	{
	return origCmdName;
	}
	
/*---------------------------------------------------------------------------------------------------
|	Used by the NxsCommandManager to control whether or not the command can execute it's controlled
|	function
*/
inline void NxsCommand::SetAllowCommandExecution(
  bool canExecute) /* true if the  controlled should be executed if the command is parsed */
	{
	canExecuteHandler = canExecute;
	}
	
/*---------------------------------------------------------------------------------------------------
|	base class version does nothing this function is called when there are non-fatal warnings to 
|	report.  Hook for NxsAutomaticallyParsedCommand to overload
*/
inline	void NxsCommand::PrintWarnings() const 
	{ 
	} 

/*---------------------------------------------------------------------------------------------------
|	Called by the NxsCommandManager to set the minimal abbreviation
*/
inline void NxsCommand::SetAbbreviation(
  const std::string & n)
	{
	abbrevCmdName = n;
	}	

/*---------------------------------------------------------------------------------------------------
|	makes sure the command can be read multiple times or hasn't been read in this block
*/
inline bool NxsCommand::CanBeRead() const
	{
	return (numTimesCmdExpected == NxsCommand::kNoLimit || numTimesRead == 0 || numTimesCmdExpected == NxsCommand::kAtLeastOnce);
	}

/*---------------------------------------------------------------------------------------------------
|	makes sure the command was read at least once in this block if it is required.
*/
inline bool NxsCommand::WasFoundCorrectNumberOfTimes() const 
	{
	if (numTimesCmdExpected == NxsCommand::kNoLimit || numTimesRead == 1)
		return true;
	if (numTimesRead == 0)
		return (numTimesCmdExpected == NxsCommand::kOnceMax);
	return (numTimesCmdExpected == NxsCommand::kAtLeastOnce);
	}

/*---------------------------------------------------------------------------------------------------
|	Determines whether or not abbreviations are allowed instead of full keyword names
*/
inline void NxsCommand::SetAllowAbbreviatedKeywords(
  bool letAbbreviate)
  	{
  	allowAbbrevKeywords = letAbbreviate;
  	}
	
inline void NxsCommand::SetDescription(const std::string & d)
	{
	descrip = d;
	}
	
inline const std::string & NxsCommand::GetDescription() const
	{
	return descrip;
	}

inline void NxsCommand::ExpectEqualsSign()
	{
	AddUnnamedSetting(NxsCommand::GetEqualsSignRequiredTokenShPtr(), false);
	}

	
#endif
