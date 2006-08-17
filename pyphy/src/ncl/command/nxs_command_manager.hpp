#ifndef NCL_NXS_COMMAND_MANAGER_H
#define NCL_NXS_COMMAND_MANAGER_H
#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/command/nxs_command_decl.hpp"
#include "pyphy/src/ncl/misc/string_extensions.hpp"
#include <deque>
class NxsToken; 
class NxsBasicListManager;
class XCommandNotFound 
	{};	//thrown by GetCommand
	
/*--------------------------------------------------------------------------------------------------------------------------
| 	Class that coordinates the parsing and execution of a single command line.
|	Used as the base class for CmdMgrNexusBlocks to provide most of the functionality for read command.
|	ProcessCommandLine() reads the input stream into tokens until a semicolon is reached, then decides which Command is 
|	being called and passes a container of tokens to that Command to parse. 
|
|	A NxsCommandManager can be told to supply help.  This means that it will check for requests for help or ? as the first 
|	token.  If either is found ProcessHelp will be called (it either lists all commands or provides help for a specific
|	command if "help commandName" was entered).  If the "supply help" option is in effect, Commands are told to skip
|	help request tokens in the middle of command strings.  If one is found the command will then display its settings
|	Note that any settings will be read if there is a ? in the middle of a command (e.g. dimensions ? ntax = 5; would 
|	set ntax to 5 and display the help)
|
|	A NxsCommandManager can also be told to provide a set command.  This allows users to specify options to a command 
|	without executing the command.  Note that if set is added, there can still be a set command that controls some options.
|	the set command provided by the NxsCommandManager will check if the second word in the command line, is a command name.
|	if it is that command will be read but not executed.  If the second word is not a command name the commandManager will
|	check for a programmers-supplied set command.
|	
|	Correct supplying of the set and help commands depends on the commands following the protocol of having all parsing
|	functionality in one function Command::ParseCommand() (this function should skip help requests when told to check
|	for them) and all execution functionality in a separate function (accessed through 
|	Command::ExecuteControlledFunction()).
|
|	Note that Commands are added as boost::shared_ptr.  It is nice to be able to use some commands outside of the 
|	scope NxsCommandManager, but for other commands letting the NxsCommandManager "own" them is ideal. 
*/
class NxsCommandManager
	{
	public:
		typedef std::pair<std::string, NxsBasicListManager *> IndexManagerInfo;
		typedef std::vector<IndexManagerInfo> VecIndexManagers;
		
		NxsCommandManager(bool abbrev = false, bool giveHelp = false, bool giveSet = false, bool giveQuery = false);
		virtual ~NxsCommandManager();
		
		//	Accessors
		//
		const VecNxsCommandShPtr & GetAllCommands() const
			{
			return commands;
			}
		const VecIndexManagers & GetCommandEnvironments() const
			{
			return indexManagers;
			}
		NxsCommandShPtr 	GetCommand(const std::string & n) const;
		bool				GetSupplyingHelpCmd() const;
		bool				GetSupplyingSetCmd() const;
		//	Modifiers
		//
		void AddCommandEnvironment(const IndexManagerInfo & mgr)
			{
			indexManagers.push_back(mgr);
			}
		NxsCommandShPtr 	AddCommand(NxsCommandShPtr);
		bool				FinishedAddingCommands();
		bool				SetAllowAbbreviations(bool letAbbrev = true, bool passOntoCommands = true);
		void				SetPrintStateToHiddenStream(bool printState = true);
		void				SetSupplyingHelpCmd(bool doAddHelp = true);
		void				SetSupplyingSetCmd(bool doAddSet = true);
		//	Utilities
		//
		void				ProcessCommandLine(NxsToken& token);
		void				ListCommandNames();
		void				ProcessHelp(NxsToken & token);
		void				ProcessAvailable(NxsToken & token);

	protected :
		typedef std::map<std::string, std::string, NxsStringNoCaseLess> AliasMap;
		typedef AliasMap * AliasMapPtr;
		bool				DetermineAbbreviations();
		NxsCommandShPtr		GetCommand(unsigned which_cmd) const;
#		if ALLOWING_MULTI_WORD_CMD_NAME
			unsigned			GetCommandIndex(NxsToken &) const;
			bool				allowMultiWordCmdNames;
#		endif
		AliasMapPtr			GetAliasMap();
		unsigned			GetCommandIndex(const std::string &n) const;
		std::string			GetCommandName(unsigned which_cmd) const;
		unsigned 			GetNumCommands() const;
		void				FinishedReadingBlock() const;
		void 				PrepareToReadNewBlock();
		bool 				PreProcessCommandLine(NxsToken &token, bool *helpReq);
		void				PrintStateToHiddenStream() const;
		void 				ProcessAliasedCommandLine(const std::string &, NxsToken &token);
		
		VecNxsCommandShPtr		commands;			/* a vector of pointers to the commands (with the index given by OakKernel's command enumeration)*/
		bool					addHelp;			/* whether or not help command should be added  */
		bool					addSet;				/* whether or not the set command should be added */
		bool					addAvailable;		/* whether or not the 'available' command should be added */
		bool					allowAbbreviations; /* whether or not commands names can be abbreviated */
		bool					readyToParse;		/* true if the manager can read command (FinishedAddingCommands() has been called with no errors ) */
		bool					printStateToHiddenStream; /// true to support printing a text version of the command state to the hidden_query stream
		NxsCommandShPtr			explicitHelp;
		std::deque<std::string> commandHistory;
		AliasMap				aliases;
		VecIndexManagers		indexManagers;
	private :
		void				Clear();
		// Copy constructor and = operator should not be used because of aliasing
		NxsCommandManager(const NxsCommandManager &);	
		NxsCommandManager& operator=(const NxsCommandManager &);
		
		friend class NxsCommandManagerBlock;
	};

inline std::map<std::string, std::string, NxsStringNoCaseLess> * NxsCommandManager::GetAliasMap()
	{
	return &aliases;
	}
	
inline unsigned NxsCommandManager::GetNumCommands() const
	{
	return (unsigned)commands.size();
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
| 	returns a pointer to the command with the name n, throws XCommandNotFound() if no command of that name is found
*/
inline NxsCommandShPtr NxsCommandManager::GetCommand(
  const std::string &n) const
  	{
  	return GetCommand(GetCommandIndex(n));
  	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	clears the list of recognized commands and flags the NxsCommandManager as being unable to read commands
*/
inline void NxsCommandManager::Clear()
	{
	readyToParse = false;
	commands.clear();
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns true if the help command is being supplies by the NxsCommandManager
*/
inline bool NxsCommandManager::GetSupplyingHelpCmd() const
	{
	return addHelp;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns true if the set command is being supplies by the NxsCommandManager
*/
inline bool NxsCommandManager::GetSupplyingSetCmd() const
	{
	return addSet;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Tells the NxsCommandManager whether or not to supply the help command
*/
inline void NxsCommandManager::SetSupplyingHelpCmd(
  bool doAddHelp)/* true to supply the help command, false to disable this feature */
	{
	addHelp = doAddHelp;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Tells the NxsCommandManager whether or not to supply the set command
*/
inline void NxsCommandManager::SetSupplyingSetCmd(
  bool doAddSet)	/* true to supply the set command, false to disable this feature */
	{
	addSet = doAddSet;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	returns a pointer to the command with index which_cmd
*/
inline NxsCommandShPtr NxsCommandManager::GetCommand(
  unsigned which_cmd) const /*index of the command from the Kernel's enum of commands) */
	{
	if (which_cmd >= commands.size())
		{
		NXS_ASSERT(0);
		throw XCommandNotFound();
		}
	return commands[(VecNxsCommandShPtr::size_type) which_cmd];
	}

inline void NxsCommandManager::SetPrintStateToHiddenStream(bool printState)
	{
	printStateToHiddenStream = printState;
	}

#endif
