#ifndef PHO_FLOOR_H
#define PHO_FLOOR_H

#include <stack>
#include "ncl/nxs_block.hpp"
#include "ncl/nxs_reader.hpp"
#include "ncl/output/nxs_output.hpp"
#include "ncl/output/nxs_table.hpp"
#include "ncl/misc/nxs_file_path.hpp"
#include "ncl/taxa/nxs_taxa_manager.hpp"	
#include "ncl/taxa/nxs_taxa_listener.hpp"
#include "ncl/trees/nxs_tree_listener.hpp"
#include "ncl/characters/nxs_char_listener.hpp"
#include "phycas/misc/long_operation_manager.hpp"
class AliasSettings;
class ExecuteSettings;
class SetSettings;
class PhoMCMC;
class DistributionManager;
class PhoSampler;

#if defined(SCM_MODULE) //POL-23Nov2004
	class SCM;
#endif

#if defined(GG_MODULE) //POL-14April2005
	class GG;
#endif

#if defined(PHYLODIVERSITY_MODULE) 
	class Phylodiversity;
#endif

#if defined(SRQ_MODULE) //POL 9-September-2004
	class SRQ;
#endif
#if defined(SUPPORT_GETTREES)
	class GetTreesSettings;
	class ClearTreesSettings;
#endif 

class PWD;
class PhoTaxaManager;
class NxsTaxaManager;
class NxsTreesManager;
class PhoTreesManager;
class PhoCharactersManager;
#if ! (NEW_NXS_BLOCK_AND_READER)
#   error NEW_NXS_BLOCK_AND_READER must be defined
#endif
/*----------------------------------------------------------------------------------------------------------------------
|	PhoFloor provides a template for creating a program that reads NEXUS data files and provides a basic command 
|	line. After compiling PhoFloor, you will already have a program that understands the following commands, either 
|	typed in at the console or provided in a PhoFloor block in a NEXUS data file (exception is the execute command,
|	which can only be entered at the console). Keywords in the descriptions below are given in uppercase, however the
|	commands themselves are case-insensitive. Lower-case indicates a parameter supplied by the user (e.g., "filename" 
|	would be replaced by the actual name of the file). Square brackets indicate optional keywords or subcommands.
|>
|	EXECUTE filename;
|	
|	LOG [options];
|	
|	  Option         Action
|	  ------------------------------------------------------
|	  FILE=filename  specifies name of log file to start
|	  START          indicates logging is to be started
|	  STOP           indicates logging is to be stopped
|	  APPEND         append to log file if it already exists
|	  REPLACE        replace log file without asking
|	
|	QUIT;
|>
|	See the Read function for details and to add other commands.
|	
|	To change the name of the program (which is also the prompt name and the name of the program's private NEXUS 
|	block), replace all occurrences of PhoFloor with the name of your program (also search for the string 
|	"PhoFloor" and replace with an appropriate string at each occurrence).
|	
|	This class handles reading and storage for the NxsReader block PhoFloor. It also serves as the main class for 
|	the program PhoFloor, acting as both a NxsReader object (in order to be capable of parsing data files) as well 
|	as a NxsBlock object (in order to be able to process commands in a PhoFloor block). 
|
|	Acting as a NxsBlock, it overrides the member functions Read and Reset, which are virtual functions in the base 
|	class NxsBlock. Acting as a NxsReader object, it overrides the member functions EnteringBlock, SkippingBlock, and 
|	NexusError.
|	
|	Adding a new data member? Don't forget to:
|~
|	o Describe it in the class header comment at the top of "PhoFloor.h"
|	o Initialize it (unless it is self-initializing) in the constructor and reinitialize it in the Reset function
|	o Describe the initial state in the constructor documentation
|	o Delete memory allocated to it in both the destructor and Reset function
|	o Report it in some way in the Report function
|~
*/
class PhoFloor
  :
# if defined(READ_PHYCAS_BLOCK)
	public NxsCommandManagerBlock, 
# endif
  public NxsReader,
  public NxsCharListener,
  public NxsTreeListener,
  public NxsTaxaListener
	{
	public:
							PhoFloor(); 
		virtual				~PhoFloor();
		
		
		void				Run(char *infile_name);
		// Nexus reading
		// 
#		if defined(READ_PHYCAS_BLOCK)
			void 				AddRecognizedCommands();
#		endif
		void 				NexusError(const std::string &msg, file_pos pos, unsigned line, unsigned col, CmdResult resCode, NxsBlock* b = NULL);
		CmdResult 			EndEncountered();
		bool 				EnteringBlock(const std::string &blockName );
		void 				ExitingBlock(const std::string &blockName );
		void 				SkippingDisabledBlock(const std::string &blockName );
		void 				SkippingBlock(const std::string &blockName );
		
		void				TaxaChanged(BaseTaxaManager *, NxsTaxaListener::TaxaChangeType);
		void				TreesChanged(NxsTreesManager *, NxsTreeListener::TreeChangeType);
		void				CharsChanged(NxsCharactersManager *, NxsCharListener::CharChangeType);

		// Access 
		//
		bool				ReadTokenStream(NxsToken &token);
		bool GetIsExiting() const
			{
			return quitNow;
			}
		void				Prompt() ;

	protected:
		bool				HandleNextCommand(const char *);
		void				PreprocessNextCommand(std::string *);
		virtual void 		PhorestYield();
		void				ResetCmdMgrNxsBlock() {}	//nothing to be done in Reset
		
		VecString			GetCommandNames();
		
#		if defined(READ_PHYCAS_BLOCK)
			CmdResult			HandleAlias(AliasSettings *);
			CmdResult			HandleCopying();
			CmdResult			HandleExecute(ExecuteSettings *);
#			if defined(SUPPORT_GETTREES)
				CmdResult			HandleGetTrees(GetTreesSettings *);
#			endif 
			CmdResult			HandleQuit();
			CmdResult			HandleSet(SetSettings *);
			
			bool 				ParseHelp(NxsToken &token);
			void 				SetupAliasCommand(NxsCommandManager *);
			void 				SetupCopyingCommand(NxsCommandManager *);
			void 				SetupExecuteCommand(NxsCommandManager *);
#			if defined(SUPPORT_GETTREES)
				void			SetupGetTrees(NxsCommandManager *);
#			endif 
			void 				SetupHelpCommand(NxsCommandManager *);
			void 				SetupQuitCommand(NxsCommandManager *);
			void 				SetupSetCommand(NxsCommandManager *);
#		endif
		std::string			next_command;	/* string that holds the most recent command */
		bool				quitNow;
		
#		if defined(SCM_MODULE) //POL-27Jan2005
			boost::shared_ptr<SCM>	scm;
#		endif

#		if defined(GG_MODULE) //POL-14April2005
			boost::shared_ptr<GG>	gg;
#		endif

	private:
		std::stack<NxsInFilePath>	fileStack;
		std::stack<NxsOutputOperationStatusID> blockOpIDStack;
		
		NxsOutputManager	  & outputMgrRef; // alias
		bool					notifyOfUnknownBlocks; /* when true a "skipping unknown XXX block" message is displayed when reading a file with unsupported blocks */
		boost::shared_ptr<PhoTaxaManager> phoTaxaMgr;
		boost::shared_ptr<PhoTreesManager> phoTreesMgr;
		boost::shared_ptr<PhoCharactersManager> phoCharactersMgr;
	}; 

/*----------------------------------------------------------------------------------------------------------------------
|	Override this virtual function in GUI applications to yield control to the operating system for purposes of
|	processing messages being directed to other applications.
*/
inline void PhoFloor::PhorestYield()
	{
	}

#endif

