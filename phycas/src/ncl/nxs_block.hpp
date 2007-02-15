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

#ifndef NCL_NXSBLOCK_H
#define NCL_NXSBLOCK_H

#include "phycas/src/ncl/nxs_defs.hpp"

#undef NCL_NXSBLOCK_SILENTLY_SKIPS_CMDS


// April 2003 changes:
//	NxsBlock is now abstract Read must be overridden.
//	NxsCommandManagerBlock is derived from NxsBlock.  It overrides Read(), but AddRecognizedCommands() now must be defined by 
//		each class derived from it.
//	To allow one NxsBlock to parse multiple types of blocks in a NEXUS fiel, NxsBlock::GetID() is virtual and is now called 
//		after the block was read.  The virtual CanReadBlockType( const std::string &) is now used to determine whether a NxsBlock object can
//		read a block in the file.  This change was mainly made so that the object that parses CHARACTERS blocks can also parse DATA blocks
//		The base class behavior of CanReadBlockType is simply to check the ID, so NxsBlocks that only read one type of NEXUS blocks shouldn't
//		need to change.
//	Report now takes the typedef NxsOutputStream instead of ostream.  see notes on nxsoutput.h for specifying output styles
//	All blocks now have space for a title
//

#include "phycas/src/ncl/command/nxs_command_manager.hpp"

class NxsReader;

/*----------------------------------------------------------------------------------------------------------------------
|	interface that encapsulates the ability to read a block of nexus commands.  
|	Pointers to NexusBlocks are added to the NxsReader object which is in charge of reading streams of nexus tokens 
|	(typically files).
|	when the NxsReader object encounters a name of a block it calls each NxsBlock's CanReadBlockType() function with the
|	block name as an argument.
|	If the object returns true, then IsEnabled() is called. 
|	If true is returned then NxsBlock::Reset(), NxsBlock::Read(), and NxsBlock::Report() are called.
|	Note that Read and Report are pure virtual functions
|	The base class keeps a pointer to NxsReader object, the block's id, and info on whether the block has been read, is enabled
|	and the strictness of the NEXUS parsing.
*/
class NxsBlock
	{
	public:
		NxsBlock(const std::string &s);
		virtual ~NxsBlock();
		
		//	Accessors/Queries
		//
		virtual bool		CanReadBlockType( const std::string &bt);
		std::string  		GetErrorMessage() const;
		virtual std::string GetID() const;
		std::string			GetTitle() const;
		bool 				IsEmpty() const;
		bool 				IsEnabled() const;
		bool 				IsUserSupplied() const;
		bool				ExitRequested() const;
		//	Modifiers
		//
		void 				Disable();
		void 				Enable();
		void				Reset();
		void 				SetNexus( NxsReader* nxsptr );
		void				SetID(const std::string &s);
		void				SetTitle(const std::string &t);

		void				SkippingCommand(const std::string &) const;

		//	Utilities
		//
		virtual bool		Read( NxsToken& token ) = 0;	
		virtual void 		Report( NxsOutputStream& out ) const;
		
	protected:
		void				ResetNxsBlockBase();
		virtual void		ResetNxsBlock() = 0;

		mutable std::string 	errorMsg;
		bool 		isEmpty;
		bool 		isEnabled;
		bool		isUserSupplied;
		bool		skipUnknownCommands;
		bool		stopReading;		/* can be set to true if a leave or quit command is encountered (if this command manager recognizes such commands) e.g. checked in NexusBlock::UseCmdMgrToRead*/
		NxsReader  *nexus;
		std::string 	blockID;		/* the type of the block */
		std::string	title;	/*individual block's name */
	};
	
#if defined (NCL_SUPPORT_OLD_NAMES)
	typedef NxsBlock NexusBlock;
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	returns the stopReading field (should be set to true if a quit, exit, or leave command is recognized by the block)
*/
inline bool NxsBlock::ExitRequested() const
	{
	return stopReading;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Uses the NxsCommandManager interface to provide the Read() functionality of NxsBlock to make it easier to write and
|	maintain new blocks.  This alleviates the need to write lots of tedious and error prone checking of user supplied 
|	options, and implements a standardized method of checking for and reporting errors to the user
|	Note derived classes must still override NxsBlock::Report and now AddRecognizedCommands() must be overridden.
|	InitializeRecognizedCommands should be called in the constructor and/or Reset of a derived class.
*/
class NxsCommandManagerBlock : public NxsCommandManager, public NxsBlock
	{
	public:
		NxsCommandManagerBlock(const std::string &s, bool abbrev = false, bool giveHelp = false, bool giveSet = false, bool giveAvailable = false);
		~NxsCommandManagerBlock();
		virtual bool CanReadBlockType(const std::string &s) 
			{
			return NxsBlock::CanReadBlockType(s);
			}
		bool				Read( NxsToken& token );	
		
	protected:
		virtual void 		AddRecognizedCommands() = 0;
		virtual CmdResult	EndEncountered();
		void 				InitializeRecognizedCommands();
		virtual void		ResetCmdMgrNxsBlock() = 0;
		void				ResetNxsBlock();
		void				SkipOneCommandLine(NxsToken &);
		virtual void 		SkippingCommand(const std::string &commandName );
		
	};

inline CmdResult NxsCommandManagerBlock::EndEncountered()
	{
	return kCmdSucceeded;
	}	

/*----------------------------------------------------------------------------------------------------------------------
|	returns true if this block can read the commands following BEGIN blockName.
|	base class version merely compares s to the idstring.
|	Can be overriden if the same NxsBlock object can read more than one type of block (e.g. Characters and Data)
 */
inline bool NxsBlock::CanReadBlockType(
  const std::string &blockName)
  	{
  	return EqualsCaseInsensitive(blockName, GetID());
  	}

/*----------------------------------------------------------------------------------------------------------------------
|	should be called in derived class's constructor and/or Reset functions if the NxsCommandManager interface is to be used
|	to read the blocks.
|	This function clears old commands and calls NxsBlock::AddRecognizedCommands() (which should be overridden to 
|	create and add the commands that apply to this block) and then NxsCommandManager::FinishedAddingCommands().
 */
inline void NxsCommandManagerBlock::InitializeRecognizedCommands()	
	{
	NxsCommandManager::Clear();
	AddRecognizedCommands();
#	if defined(NDEBUG)
		FinishedAddingCommands();
#	else
		bool check = FinishedAddingCommands();
		NXS_ASSERT(check);
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Does nothing.
 */
inline NxsCommandManagerBlock::~NxsCommandManagerBlock()
	{
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the title std::string.
 */
inline std::string NxsBlock::GetTitle() const
	{
	return title;
	}
	
inline void NxsCommandManagerBlock::SkippingCommand(
  const std::string &unknownCmd) 	/* the nexus token that was not recognized as a command */
	{
	NxsBlock::SkippingCommand(unknownCmd);
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Calls NxsBlock::Reset() and NxsCommandManager::PrepareToReadNewBlock()
*/
inline void NxsCommandManagerBlock::ResetNxsBlock()
	{
	PrepareToReadNewBlock();
	ResetCmdMgrNxsBlock();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	sets the id string (which, by default, is used in the CanReadBlockType() function to determine whether or not this
|	block can read a group of command following BEGIN XXX)
*/
inline void NxsBlock::SetID(
  const std::string &s)	/* the block ID, for example "TAXA" */
	{
	blockID = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	sets the title (name) of this individual block
*/
inline void NxsBlock::SetTitle(
  const std::string &s)	/* the block title, for example "FirstTAXABlock" */
	{
	title = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This virtual function provides a hook for displaying a brief report of the contents of the block.
*/
inline void NxsBlock::Report(
  NxsOutputStream &) const /* the output stream to which the report is sent */
	{
	}
		

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the nexus data member of the NxsBlock object to 'nxsptr'.
*/
inline void NxsBlock::SetNexus(
  NxsReader *nxsptr)	/* pointer to a NxsReader object */
	{
	nexus = nxsptr;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if Read function has not been called since the last Reset. This base class version simply returns the 
|	value of the data member isEmpty. If you derive a new block class from NxsBlock, be sure to set isEmpty to true in 
|	your Reset function and isEmpty to false in your Read function.
*/
inline bool NxsBlock::IsEmpty() const
	{
	return isEmpty;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Should be overridden but still call the base class version.
| 	Reset is Called by the NxsReader object just prior to calling the block object's Read function.
|	 
*/
inline void NxsBlock::Reset()
	{
	ResetNxsBlockBase();
	ResetNxsBlock();
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	NxsBlock function resets the errorMsg to an empty string and flags the block  as empty
 */
inline void NxsBlock::ResetNxsBlockBase()
	{
	errorMsg.clear();
	isEmpty = true;
	isUserSupplied = false;
	stopReading = false;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Does nothing.
 */
inline NxsBlock::~NxsBlock()
	{
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of isEnabled to false. A NxsBlock can be disabled (by calling this method) if blocks of that type
|	are to be skipped during execution of the NEXUS file. If a disabled block is encountered, the virtual
|	NxsReader::SkippingDisabledBlock function is called, giving your application the opportunity to inform the user
|	that a block was skipped.
*/
inline void NxsBlock::Disable()
	{
	isEnabled = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of isEnabled to true. A NxsBlock can be disabled (by calling Disable) if blocks of that type are to
|	be skipped during execution of the NEXUS file. If a disabled block is encountered, the virtual 
|	NxsReader::SkippingDisabledBlock function is called, giving your application the opportunity to inform the user
|	that a block was skipped.
*/
inline void NxsBlock::Enable()
	{
	isEnabled = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of isEnabled, which can be controlled through use of the Enable and Disable member functions. A 
|	NxsBlock should be disabled if blocks of that type are to be skipped during execution of the NEXUS file. If a 
|	disabled block is encountered, the virtual NxsReader::SkippingDisabledBlock function is called, giving your 
|	application the opportunity to inform the user that a block was skipped.
*/
inline bool NxsBlock::IsEnabled() const
	{
	return isEnabled;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of isUserSupplied, which is true if and only if this block's Read function is called to process a 
|	block of this type appearing in a data file. This is useful because in some cases, a block object may be created 
|	internally (e.g. a NxsTaxaBlock may be populated using taxon names provided in a DATA block), and such blocks do 
|	not require permission from the user to delete data stored therein.
*/
inline bool NxsBlock::IsUserSupplied() const
	{
	return isUserSupplied;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the id std::string.
 */
inline std::string NxsBlock::GetID() const
	{
	return blockID;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns a copy of the error message that is associated the NxsBlock being unable to read a block.
|	if an NxsException exception is thrown with no message, this function is used to see if the NxsBlock::errorMsg indicates
|	what type of error occurred.
*/
inline std::string NxsBlock::GetErrorMessage() const
	{
	return errorMsg;
	}

#endif


