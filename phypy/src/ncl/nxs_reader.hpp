#ifndef NCL_NXSREADER_H
#define NCL_NXSREADER_H
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/output/nxs_output.hpp"
#if NEW_NXS_BLOCK_AND_READER

#include <boost/weak_ptr.hpp>
#include <list>
class NxsBlock;
typedef boost::weak_ptr<NxsBlock> NxsBlockWeakPtr;
typedef boost::shared_ptr<NxsBlock> NxsBlockPtrShared;

class NxsToken;
//	changes
//	blockList is now a list of pointers as opposed to simple pointer (the next pointer in the NxsBlock is removed)
//	Execute interacts with new NexusTokenizer
//		To determine whether or not a block can read a group or commands.  Now calls NxsBlock::CanReadBlockType(std::string) instead of comparing to the idStr (more flexible)
//	const-correct
//	Now accepts smart pointers to blocks (NxsBlockWeakPtr), if user doesn't want to be in charge of deleting blocks.
//	interface changes
//  deprecated
//		PositionInBlockList (unnecessary)
//		Reassign  (use Detach then Add instead of Reassign)
//	allowMissingInEquate and allowPunctuationInNames not supported (were they ever?)
//
/*----------------------------------------------------------------------------------------------------------------------
|	Class that coordinates the reading of a stream of Nexus tokens (through the Execute() function).  
|	NxsReader stores a list of the NexusBlocks objects that can read the file  (pointers to these objects are added and 
|	removed using the Add and Detach functions).
|	Several virtual functions called in the Execute function allow flexibility in user-interface and pre- and post-
|	processing of the blocks.
|	
|	To allow derived classes to prepare for or absorb the input of more data:
|	ExecuteStarting and ExecuteStopping are called before and after the stream of tokens is read.
|	EnteringBlock and ExitingBlock are called before and after each block (with the name of the block as a parameter).
|	
|	To allow derived classes to interact with the user:
|	SkippingDisabledBlock is called when a block is disabled but its name is encountered in the file.  
|	SkippingBlock is called when a block name in the file is does not match any block that has been added.
|	OutputComment() is called when an [! ... ] style comment is encountered
|	NexusError() is called when any error is found in the stream of tokens (with arguments describing the
|	error and its location).
|
|	Note: all blocks added are treated as aliases.  The NxsReader object does not delete any of them in its destructor
*/
class NxsReader
	{

	public:
		/*  error codes for UseNexusBlockToRead() */
		enum nxs_block_error_codes
			{	block_read, 		/*  no error*/
				disabled_block, 	/* NxsBlock::IsEnabled returned false (SkippingDisabledBlock already called)*/
				could_not_enter, 	/* NxsReader::EnteringBlock() returned false*/
				early_eof, 			/*  unexpected end of file in the block (NexusError has already been called to report the error)*/
				error_in_block,		/* NxsException error thrown in reading the block (NexusError has already been called to report the error) */
				block_not_found,		/* the default setting of the errCode will be overwritten if any attempt is made to read the block */
				quit_or_leave
			}; 	

		virtual ~NxsReader();

		inline STATIC_CONST const char   *NCLNameAndVersion();
		inline STATIC_CONST const char   *NCLCopyrightNotice();
		inline STATIC_CONST const char   *NCLHomePageURL();
		

		// Accessors
		//
				bool	BlockListEmpty() const;
				
		// Modifiers
		//
#		if defined (NCL_USING_BLOCK_SUPPLIERS)
			void 			AddBlockSupplier(NxsBlockSupplierShPtr newBlockMgr );
#		endif
		void 			Add(NxsBlock *newBlock);
		void 			Add(NxsBlockWeakPtr  newBlock);
		void			ResetAllBlocks();
		
		void 			Detach(NxsBlock* block );
		void 			Detach(NxsBlockWeakPtr  block );
		
		// Utilities
		//
		std::vector<bool> 		DisableAllBlocksExcept(const std::string &activeBlockID);
		std::vector<bool> 		DisableAllBlocksExcept(const VecString &activeBlockIDs);
		bool 					Execute( NxsToken & token, bool notifyStartStop = true );
		nxs_block_error_codes 	UseNexusBlockToRead(NxsBlock *currBlock, NxsToken &token);
		void					RestoreEnabledStatus(const std::vector<bool> &); //assumes that the # and order of block readers hasn't changed

		// Debugging
		//
		virtual void 	DebugReportBlock( const NxsBlock& nexusBlock ) const;
		
		// Polymorphic Functions
		//
		virtual bool 	EnteringBlock(const std::string &blockName );
		virtual void 	ExecuteStarting();
		virtual void 	ExecuteStopping();
		virtual void 	ExitingBlock(const std::string &blockName );
		virtual void 	NexusError(const std::string &msg, file_pos pos, unsigned line, unsigned col, CmdResult resCode, NxsBlock* b = NULL);
		virtual void 	SkippingDisabledBlock(const std::string &blockName );
		virtual void 	SkippingBlock(const std::string &blockName );

		
	protected:
		typedef std::list<NxsBlock*> 			ListOfBlockPtrs;
		typedef std::list<NxsBlockWeakPtr> 		ListOfWeakBlockPtrs;
		ListOfBlockPtrs 		blockPtrs;				/* set of all bare pointers to block readers that can be used to read the stream of tokens */
		ListOfWeakBlockPtrs 	weakBlockPtrs;		/* set of all smart pointers to block readers that can be used to read the stream of tokens */
#		if defined(NCL_USING_BLOCK_SUPPLIERS)
			typedef list<NxsBlockSupplierShPtr> listOfBlockMgrs;
			listOfBlockMgrs 		blockMgrs;				/* set of all objects that can supply block readers that can be used to read the stream of tokens */
#		endif		
	};

#if defined (NCL_SUPPORT_OLD_NAMES)
	typedef NxsReader Nexus;
#endif

/*----------------------------------------------------------------------------------------------------------------------
| 	Does nothing.  NOT responsible for deleting the NxsBlock objects
*/
inline NxsReader::~NxsReader()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	returns true if NxsBlock objects have been added (or all have been removed using the Detach function).
 */
inline bool NxsReader::BlockListEmpty() const
	{
	return (blockPtrs.empty() && weakBlockPtrs.empty());
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	This function was created for purposes of debugging a new NxsBlock.
| 	This version does nothing; to create an active DebugReportBlock function,
| 	override this version in the derived class and call the Report function
| 	of nexusBlock.  This function is called whenever the main NxsReader Execute function
| 	encounters the [&spillall] command comment between blocks in the data file.
| 	The Execute function goes through all blocks and passes them, in turn, to this
| 	DebugReportBlock function so that their contents are displayed. Placing the
| 	[&spillall] command comment between different versions of a block allows
| 	multiple blocks of the same type to be tested using one long data file.
| 	Say you are interested in testing whether the normal, transpose, and
| 	interleave format of a matrix can all be read correctly.  If you put three
| 	versions of the block in the data file one after the other, the second one
| 	will wipe out the first, and the third one will wipe out the second, unless
| 	you have a way to report on each one before the next one is read.  This
| 	function provides that ability.
*/
inline void NxsReader::DebugReportBlock(const NxsBlock& /* nexusBlock */)const
{
   // Derive me and uncomment the following line in the derived version
   // nexusBlock.Report(out);
	// Note that your derived NxsReader object must have its own ostream object (out)
}

/*----------------------------------------------------------------------------------------------------------------------
|	Called by the NxsReader object when a block named `blockName' is entered. Allows derived class overriding this
|	function to notify user of progress in parsing the NEXUS file. Also gives program the opportunity to ask user if it
|	is ok to purge data currently contained in this block. If user is asked whether existing data should be deleted, and
|	the answer comes back no, then then the overrided function should return false, otherwise it should return true.
|	This (base class) version always returns true.
*/
inline bool NxsReader::EnteringBlock(
  const std::string &blockName)	/* the name of the block just entered */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(blockName)
#	endif
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called by the NxsReader object when a block named `blockName' is being exited. Allows derived class overriding this
|	function to notify user of progress in parsing the NEXUS file.
*/
inline void NxsReader::ExitingBlock(
  const std::string &blockName)	/* the name of the block being exited */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(blockName)
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when an error is encountered in a NEXUS file. Allows program to give user details of the error as well as 
|	the precise location of the error.
*/
inline void NxsReader::NexusError(
  const std::string &msg,	/* the error message to be displayed */
  file_pos	pos,	/* the current file position */
  unsigned	line,	/* the current file line */
  unsigned	col,
  CmdResult ,
  NxsBlock *) 		/* the current column within the current file line */
	{
#	if defined (NCL_USE_STD_OUTPUT)
		NxsOutputManager localOM; //@ we need to add a NxsOutputManager to the reader
		*localOM.GetErrorStreamPtr() << msg << "(at line = " << line << ", col = " << col << ')' << std::endl;
#		if defined(HAVE_PRAGMA_UNUSED)
#			pragma unused(pos)
#		endif
#	else
#		if defined(HAVE_PRAGMA_UNUSED)
#			pragma unused(msg, pos, line, col)
#		endif
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function is called when an unknown block named `blockName' is about to be skipped. Override this pure virtual
|	function to provide an indication of progress as the NEXUS file is being read.
*/
inline void NxsReader::SkippingBlock(
  const std::string &blockName)	/* the name of the block being skipped */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(blockName)
#	endif
	
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function is called when a disabled block named `blockName' is encountered in a NEXUS data file being executed.
|	Override this pure virtual function to handle this event in an appropriate manner. For example, the program may 
|	wish to inform the user that a data block was encountered in what is supposed to be a tree file.
*/
inline void NxsReader::SkippingDisabledBlock(
  const std::string &blockName)	/* the name of the disabled block being skipped */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(blockName)
#	endif
	}




#else //NEW_NXS_BLOCK_AND_READER
/*----------------------------------------------------------------------------------------------------------------------
|	This is the class that orchestrates the reading of a NEXUS data file. An object of this class should be created, 
|	and objects of any block classes that are expected to be needed should be added to `blockList' using the Add 
|	member function. The Execute member function is then called, which reads the data file until encountering a block 
|	name, at which point the correct block is looked up in `blockList' and that object's Read method called. 
*/
class NxsReader
	{
	public:
		enum	NxsTolerateFlags	/* Flags used with data member tolerate used to allow some flexibility with respect to the NEXUS format */
			{
			allowMissingInEquate	= 0x0001,	/* if set, equate symbols are allowed for missing data symbol */
			allowPunctuationInNames	= 0x0002	/* if set, some punctuation is allowed within tokens representing labels for taxa, characters, and sets */
			};

						NxsReader();
		virtual			~NxsReader();

		bool			BlockListEmpty();
		unsigned		PositionInBlockList(NxsBlock *b);
		void			Add(NxsBlock *newBlock);
		void			Detach(NxsBlock *newBlock);
		void			Reassign(NxsBlock *oldb, NxsBlock *newb);
		void			Execute(NxsToken& token, bool notifyStartStop = true);

		virtual void	DebugReportBlock(NxsBlock &nexusBlock);

		inline const char	*NCLNameAndVersion();
		inline const char	*NCLCopyrightNotice();
		inline const char	*NCLHomePageURL();

		virtual void	ExecuteStarting();
		virtual void	ExecuteStopping();

		virtual bool	EnteringBlock(std::string blockName);
		virtual void	ExitingBlock(std::string blockName);

		virtual void	OutputComment(const std::string &comment);

		virtual void	NexusError(std::string msg, file_pos pos, unsigned line, unsigned col, CmdResult c);

		virtual void	SkippingDisabledBlock(std::string blockName);
		virtual void	SkippingBlock(std::string blockName);

	protected:

		NxsBlock		*blockList;	/* pointer to first block in list of blocks */
		NxsBlock		*currBlock;	/* pointer to current block in list of blocks */
	};

typedef NxsBlock NexusBlock;
typedef NxsReader Nexus;

/*----------------------------------------------------------------------------------------------------------------------
|	Called by the NxsReader object when a block named `blockName' is entered. Allows derived class overriding this
|	function to notify user of progress in parsing the NEXUS file. Also gives program the opportunity to ask user if it
|	is ok to purge data currently contained in this block. If user is asked whether existing data should be deleted, and
|	the answer comes back no, then then the overrided function should return false, otherwise it should return true.
|	This (base class) version always returns true.
*/
inline bool NxsReader::EnteringBlock(
  std::string blockName)	/* the name of the block just entered */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(blockName)
#	endif
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called by the NxsReader object when a block named `blockName' is being exited. Allows derived class overriding this
|	function to notify user of progress in parsing the NEXUS file.
*/
inline void NxsReader::ExitingBlock(
  std::string blockName)	/* the name of the block being exited */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(blockName)
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when an error is encountered in a NEXUS file. Allows program to give user details of the error as well as 
|	the precise location of the error.
*/
inline void NxsReader::NexusError(
  std::string msg,	/* the error message to be displayed */
  file_pos	pos,	/* the current file position */
  unsigned	line,	/* the current file line */
  unsigned	col,	/* the current column within the current file line */
  CmdResult c)
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(msg, pos, line, col)
#	endif
	}
/*----------------------------------------------------------------------------------------------------------------------
|	This function is called when an unknown block named `blockName' is about to be skipped. Override this pure virtual
|	function to provide an indication of progress as the NEXUS file is being read.
*/
inline void NxsReader::SkippingBlock(
  std::string blockName)	/* the name of the block being skipped */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(blockName)
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function is called when a disabled block named `blockName' is encountered in a NEXUS data file being executed.
|	Override this pure virtual function to handle this event in an appropriate manner. For example, the program may 
|	wish to inform the user that a data block was encountered in what is supposed to be a tree file.
*/
inline void NxsReader::SkippingDisabledBlock(
  std::string blockName)	/* the name of the disabled block being skipped */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(blockName)
#	endif
	}



#endif //NEW_NXS_BLOCK_AND_READER

/*----------------------------------------------------------------------------------------------------------------------
|	Called just after Execute member function reads the opening "#NEXUS" token in a NEXUS data file. Override this 
|	virtual base class function if your application needs to do anything at this point in the execution of a NEXUS data
|	file (e.g. good opportunity to pop up a dialog box showing progress). Be sure to call the Execute function with the
|	`notifyStartStop' argument set to true, otherwise ExecuteStarting will not be called.
|	
*/
inline void NxsReader::ExecuteStarting()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when Execute member function encounters the end of the NEXUS data file, or the special comment [&LEAVE] is
|	found between NEXUS blocks. Override this virtual base class function if your application needs to do anything at 
|	this point in the execution of a NEXUS data file (e.g. good opportunity to hide or destroy a dialog box showing 
|	progress). Be sure to call the Execute function with the `notifyStartStop' argument set to true, otherwise 
|	ExecuteStopping will not be called.
*/
inline void NxsReader::ExecuteStopping()
	{
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string containing the copyright notice for the NxsReader Class Library, useful for reporting the use of 
|	this library by programs that interact with the user.
*/
inline const char *NxsReader::NCLCopyrightNotice()
	{
	return NCL_COPYRIGHT;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string containing the URL for the NxsReader Class Library internet home page.
*/
inline const char *NxsReader::NCLHomePageURL()
	{
	return NCL_HOMEPAGEURL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string containing the name and current version of the NxsReader Class Library, useful for reporting the 
|	use of this library by programs that interact with the user.
*/
inline const char *NxsReader::NCLNameAndVersion()
	{
	return NCL_NAME_AND_VERSION;
	}


#endif

