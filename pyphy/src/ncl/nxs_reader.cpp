//#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#if NEW_NXS_BLOCK_AND_READER
#include "ncl/nxs_reader.hpp"
#include "ncl/nxs_block.hpp"
#include "ncl/nxs_token.hpp"
#include "ncl/nxs_exception.hpp"
using std::vector;
using std::string;
void DisableBlockIfItCannotReadID(NxsBlock *nxsBlockReader, const VecString &activeBlockIDs);

/*----------------------------------------------------------------------------------------------------------------------
| 	Disables any block that cannot read the activeBlockID id.  
|	returns the enabled status before the command was executed (in case you need to restore)
*/
vector<bool> NxsReader::DisableAllBlocksExcept(const string &activeBlockID)
	{
	VecString v(1, activeBlockID);
	return DisableAllBlocksExcept(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Enables blocks with a true in the vector, Disables those with false.
|	When used with DisableAllBlocksExcept this assumes that the # and order of block readers hasn't changed
*/
void NxsReader::RestoreEnabledStatus(const vector<bool> &enableStatus)
	{
	NxsBlock *nxsBlockReader;
	vector<bool>::const_iterator esIt = enableStatus.begin();
	for (ListOfBlockPtrs::iterator currIt = blockPtrs.begin(); esIt != enableStatus.end() && currIt != blockPtrs.end(); ++currIt, ++esIt)
		{
		nxsBlockReader = *currIt;
		if (*esIt)
			nxsBlockReader->Enable();
		else
			nxsBlockReader->Disable();
		}
	NxsBlockPtrShared temp = NxsBlockPtrShared();	//If we are get the pointer from a weak pointer, this shared pointer will keep the block around while it is read the block
	for(ListOfWeakBlockPtrs::iterator currIt = weakBlockPtrs.begin(); esIt != enableStatus.end() && currIt != weakBlockPtrs.end(); ++currIt, ++esIt)
		{
		temp =  boost::make_shared(*currIt); 
		nxsBlockReader = temp.get();
		if (*esIt)
			nxsBlockReader->Enable();
		else
			nxsBlockReader->Disable();
		}
	} 

void DisableBlockIfItCannotReadID(NxsBlock *nxsBlockReader, const VecString &activeBlockIDs)
	{
	bool canRead = false;
	for (VecString::const_iterator abiIt = activeBlockIDs.begin(); abiIt != activeBlockIDs.end(); ++abiIt)
		{
		if (nxsBlockReader->CanReadBlockType(*abiIt))
			{
			canRead = true;
			break;
			}
		}
	if (!canRead)
		nxsBlockReader->Disable();
	}
	
/*----------------------------------------------------------------------------------------------------------------------
| 	Disables any block that cannot read the any of the id's in activeBlockIDs.  
|	returns the enabled status before the command was executed (in case you need to restore)
*/
vector<bool> NxsReader::DisableAllBlocksExcept(const VecString &activeBlockIDs)
	{
	NxsBlock *nxsBlockReader;
	vector<bool> retVec;
	for(ListOfBlockPtrs::iterator currIt = blockPtrs.begin(); currIt != blockPtrs.end(); ++currIt)
		{
		nxsBlockReader = *currIt;
		retVec.push_back(nxsBlockReader->IsEnabled());
		DisableBlockIfItCannotReadID(nxsBlockReader, activeBlockIDs);
		}
	NxsBlockPtrShared temp = NxsBlockPtrShared();	//If we are get the pointer from a weak pointer, this shared pointer will keep the block around while it is read the block
	for(ListOfWeakBlockPtrs::iterator currIt = weakBlockPtrs.begin(); currIt != weakBlockPtrs.end(); ++currIt)
		{
		temp =  boost::make_shared(*currIt); 
		nxsBlockReader = temp.get();
		retVec.push_back(nxsBlockReader->IsEnabled());
		DisableBlockIfItCannotReadID(nxsBlockReader, activeBlockIDs);
		}
	return retVec;
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Adds newBlock to the list of NxsBlock objects that can be used to read files and calls SetNexus method of newBlock
| 	to inform newBlock of the NxsReader object that now calls it.
*/
void NxsReader::Add(
  NxsBlockWeakPtr newBlock)
	{
	assert(boost::make_shared(newBlock) != NxsBlockPtrShared());
	boost::make_shared(newBlock)->SetNexus(this);
	weakBlockPtrs.push_back(newBlock);
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Adds newBlock to the list of NxsBlock objects that can be used to read files and calls SetNexus method of newBlock
| 	to inform newBlock of the NxsReader object that now calls it.
 */
void NxsReader::Add(
  NxsBlock *newBlock)
	{
	assert(newBlock != NULL);
	newBlock->SetNexus(this);
	blockPtrs.push_back(newBlock);
	}

bool	EqualsBare(NxsBlockWeakPtr w, NxsBlock *b);
bool	EqualsBare(NxsBlockWeakPtr w, NxsBlock *b)
	{
	return (b == boost::make_shared(w).get());
	}
/*----------------------------------------------------------------------------------------------------------------------
| 	Adds newBlock to the list of NxsBlock objects that can be used to read files and calls SetNexus method of newBlock
| 	to inform newBlock of the NxsReader object that now calls it.
 */
void NxsReader::Detach(
  NxsBlockWeakPtr newBlock)
	{
	NxsBlockPtrShared p = boost::make_shared(newBlock);
	if (p != NxsBlockPtrShared())
		Detach(p.get());
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Adds newBlock to the list of NxsBlock objects that can be used to read files and calls SetNexus method of newBlock
| 	to inform newBlock of the NxsReader object that now calls it.
 */
void NxsReader::Detach(
  NxsBlock *newBlock)
	{
	assert(newBlock != NULL);
	newBlock->SetNexus(NULL);
	blockPtrs.remove(newBlock);
	//	also remove the shared block pointer to the object
	//
	for (ListOfWeakBlockPtrs::iterator sbIt = weakBlockPtrs.begin(); sbIt != weakBlockPtrs.end(); )
		{
		NxsBlockPtrShared p = boost::make_shared(*sbIt);
		if (p.get() == newBlock)
			sbIt = weakBlockPtrs.erase(sbIt);
		else
			++sbIt;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Calls the reset function for every NxsBlock object that NxsReader contains.  
*/
void NxsReader::ResetAllBlocks()
	{
	for (ListOfBlockPtrs::iterator bIt = blockPtrs.begin(); bIt != blockPtrs.end(); ++bIt)
		(*bIt)->Reset();
	for (ListOfWeakBlockPtrs::iterator sbIt = weakBlockPtrs.begin(); sbIt != weakBlockPtrs.end(); ++sbIt)
		boost::make_shared(*sbIt)->Reset();
	
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Reads the Nexus data file from the input stream provided by token.  
|	This function is responsible for reading up to the name of a each block.  Once it has read a block name, it searches
|	through the NxsBlock objects that have been Added for one that returns true for the function CanReadBlockType().
|	(if no NxsBlocks are found. the NxsBlockSuppliers are asked if they can generate a NxsBlock to read the block)
| 	If a NxsBlock reader is found NxsReader::UseNexusBlockToRead() is called with that block as an argument; 
|	if not the block is skipped (and SkippingBlock() is called).
|	This function also implements reading the command comments LEAVE and SHOWALL between blocks.
| 	The notifyStartStop argument is provided in case you do not wish the ExecuteStart and ExecuteStop functions to be 
|	called.  These functions are primarily used for creating and destroying a dialog box to show progress, and nested 
|	Execute calls can thus cause problems (e.g., a dialog box is destroyed when the inner Execute calls ExecuteStop and 
|	the outer Execute still expects the dialog box to be available).  Specifying notifyStartStop false for all the 
|	nested Execute calls thus allows the outermost Execute call to control creation and destruction of the dialog box.
|	returns true if the end of the file is reached without error.
|
|	Errors in the parsing of NxsBlocks generate NxsExceptions.  These are caught in this function, and NexusError()
|	is called to report the error, then false is returned.  Thus, this function doesn't throw any NxsExceptions.
 */
bool NxsReader::Execute(
  NxsToken& token, 
  bool notifyStartStop /* = true */)
	{
	//	Check for valid stream with the first token being #NEXUS
	//
	string errormsg;
	try 
		{
		token.ReadToken();
		}
	catch(NxsException & x) 
		{
		NexusError(x.msg, 0, 0, 0, x.cResultCode);
		return false;
		}
	if(!token.Equals("#NEXUS")) 
		{
		errormsg << "Expecting #NEXUS to be the first token in the file, but found " << token.GetTokenReference() << " instead";
		NexusError(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn(), kCmdFailedGenerateMessage);
		return false;
		}
	if(notifyStartStop)
		ExecuteStarting();
	bool leaveCommandEncountered = false;
	//	read through the file until one of the following: EOF, LEAVE command comment, or an error in a block.
	//	errors result in returning false.  All NxsExceptions are caught, funnelled through NexusError(), and result in "return false;"
	//
	for(;;)
		{
		bool cmdCommentFound = token.ReadCommandCommentOrToken();
		if(token.AtEOF())
			break;
		if (cmdCommentFound)
			{
			unsigned commandCommentIndex = token.GetNextCommentIndex('&', 0);
			if (commandCommentIndex != UINT_MAX)
				{
				NxsComment cmdComment = token.GetComment(commandCommentIndex);
				if (cmdComment.GetLocationInToken() < 1)
					{
					if (StrEquals(cmdComment.GetCommentText(), "SHOWALL", ncl::kStringNoRespectCase))
						{
						ListOfBlockPtrs::iterator currIt = blockPtrs.begin();
						for(; currIt != blockPtrs.end(); ++currIt)
							DebugReportBlock(**currIt);
						ListOfWeakBlockPtrs::iterator scurrIt = weakBlockPtrs.begin();
						NxsBlockPtrShared p = NxsBlockPtrShared();
						for(; scurrIt != weakBlockPtrs.end(); ++scurrIt)
							{
							p = boost::make_shared(*scurrIt);
							DebugReportBlock(*p);
							}
						}
					else if (StrEquals(cmdComment.GetCommentText(), "LEAVE", ncl::kStringNoRespectCase))
						{
						leaveCommandEncountered = true;
						break;
						}
					}
				}
			if (leaveCommandEncountered)
				break;
			}
		else
			{
			//	Expecting Begin BLOCK_NAME;
			//
			if (token.Equals("BEGIN"))
				{
				token.SetEOFAllowed(false);
				token.ReadToken();
				bool blockFound = false;
				NxsBlock * nxsBlockReader = NULL;
				//	Look through all of the blocks readers that have been added to the NxsReader.  If one responds true to 
				//	CanReadBlockType(), use it to read the block 
				//
				for(ListOfBlockPtrs::iterator currIt = blockPtrs.begin(); !blockFound && currIt != blockPtrs.end(); ++currIt)
					{
					nxsBlockReader = *currIt;
					if(nxsBlockReader->CanReadBlockType(token.GetTokenReference()))
						blockFound = true;
					}

				// If we get the pointer from a weak pointer, this shared pointer will keep the block around while it reads the block
				NxsBlockPtrShared temp = NxsBlockPtrShared();	
				
				if (!blockFound)
					{
					for(ListOfWeakBlockPtrs::iterator currIt = weakBlockPtrs.begin(); currIt != weakBlockPtrs.end(); ++currIt)
						{
						temp =  boost::make_shared(*currIt); 
						nxsBlockReader = temp.get();
						if(nxsBlockReader != NULL && nxsBlockReader->CanReadBlockType(token.GetTokenReference()))
							{
							blockFound = true;
							break;
							}
						}
					}
				
				int errCode = block_not_found;
				
				if (blockFound)
					{
					// A Block reader was found (or provided). Call UseNexusBlockToRead to read the info.
					// Return false for any fatal errors (disabled/unknown blocks are skipped and are 
					// not treated as fatal errors)
					//
					assert(nxsBlockReader != NULL);
					errCode = UseNexusBlockToRead(nxsBlockReader, token);
					if (errCode == quit_or_leave)
						return true;
					if (errCode == early_eof || errCode == error_in_block)
						return false;
					
					}
					
				if(errCode != block_read )
					{
					//	Unknown block or disabled block
					//
					string currBlock = token.GetTokenReference();
					try
						{
						if(errCode != disabled_block) 
							SkippingBlock(currBlock);
						for(;;)
							{
							token.ReadToken();
							if (token.Equals("END") || token.Equals("ENDBLOCK")) 
								{
								token.ReadToken();
								if(token.GetTokenReference() != ';') 
									{
									errormsg << "Expecting ';' after END or ENDBLOCK command, but found " << token.GetTokenReference() << " instead";
									NexusError(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn(), kCmdFailedGenerateMessage);
									return false;
									}
								break;
								}
							}
						}
					catch (NxsX_UnexpectedEOF & xeof)
						{
						xeof.msg << " in " << currBlock << " block.";
						NexusError(xeof.msg, xeof.pos, xeof.line, xeof.col, xeof.cResultCode, nxsBlockReader);
						return false;			
						}
					} // if token not found amongst known block IDs
				token.SetEOFAllowed(true);
				} // if token equals BEGIN
			else
				{
				errormsg << "Expecting BEGIN and then a block name, but found " << token.GetTokenReference() << " instead";
				NexusError(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn(), kCmdFailedGenerateMessage);
				return false;
				}
			}
		} 
	//	File completed (or [&LEAVE] command comment) without errors
	//
	if(notifyStartStop)
		ExecuteStopping();
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	This function is called when a NxsBlock (currBlock) is found to read the stream (token).  
|	The function calls the following functions (bailing out if any of them return an error code): 
|>
|	NxsBlock:: IsEnabled(), 
|	NxsReader::EnteringBlock(),
|	NxsBlock::Reset(),
|	NxsBlock::Read(), 
|	NxsBlock::ExitRequested(), 
|	NxsReader::ExitingBlock()
|>
|	All `NxsException's are caught and dealt with (resulting in calls to NexusError()).  
|	returns a value from the nxs_block_error_codes enum
*/
NxsReader::nxs_block_error_codes NxsReader::UseNexusBlockToRead(
  NxsBlock *currBlock,
  NxsToken &token)
	{
	if(!currBlock->IsEnabled()) 
		{
		//	hook to notify the user that a known, but disabled, block is being skipped
		//
		SkippingDisabledBlock(token.GetTokenReference());
		return disabled_block;
		}
	if(!EnteringBlock(currBlock->GetID())) 
		return could_not_enter;
	currBlock->Reset();
	try 
		{
		if (!currBlock->Read(token))
			return error_in_block;	// NOTE returning false from Read should mean that the error has been reported !!
		if (currBlock->ExitRequested())
			return quit_or_leave;
		ExitingBlock(currBlock->GetID());
		}
	catch (NxsX_UnexpectedEOF & xeof)
		{
		xeof.msg << " in " << currBlock->GetID() << " block.";
		NexusError(xeof.msg, xeof.pos, xeof.line, xeof.col, xeof.cResultCode,  currBlock);
		currBlock->Reset();
		return early_eof;
		}
	catch(NxsException & x) 
		{
		string eMsg(x.msg);
		if (eMsg.empty())
			eMsg = currBlock->GetErrorMessage();
		if ((unsigned) x.pos == 0)
			NexusError(eMsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn(), x.cResultCode, currBlock);
		else
			NexusError(eMsg, x.pos, x.line, x.col, x.cResultCode,  currBlock);
		currBlock->Reset();
		return error_in_block;
		}
	return block_read;
	}

#else //NEW_NXS_BLOCK_AND_READER
#include "ncl.h"

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes both `blockList' and `currBlock' to NULL.
*/
NxsReader::NxsReader()
	{
	blockList	= NULL;
	currBlock	= NULL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Nothing to be done.
*/
NxsReader::~NxsReader()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds `newBlock' to the end of the list of NxsBlock objects growing from `blockList'. If `blockList' points to NULL,
|	this function sets `blockList' to point to `newBlock'. Calls SetNexus method of `newBlock' to inform `newBlock' of
|	the NxsReader object that now owns it. This is useful when the `newBlock' object needs to communicate with the 
|	outside world through the NxsReader object, such as when it issues progress reports as it is reading the contents
|	of its block.
*/
void NxsReader::Add(
  NxsBlock *newBlock)	/* a pointer to an existing block object */
	{
	assert(newBlock != NULL);

	newBlock->SetNexus(this);

	if (!blockList)
		blockList = newBlock;
	else
		{
		// Add new block to end of list
		//
		NxsBlock *curr;
		for (curr = blockList; curr && curr->next;)
			curr = curr->next;
		assert(curr && !curr->next);
		curr->next = newBlock;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns position (first block has position 0) of block `b' in `blockList'. Returns UINT_MAX if `b' cannot be found
|	in `blockList'.
*/
unsigned NxsReader::PositionInBlockList(
  NxsBlock *b)	/* a pointer to an existing block object */
	{
	unsigned pos = 0;
	NxsBlock *curr = blockList;

	for (;;)
		{
		if (curr == NULL || curr == b)
			break;
		++pos;
		curr = curr->next;
		}

	if (curr == NULL)
		pos = UINT_MAX;

	return pos;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reassign should be called if a block (`oldb') is about to be deleted (perhaps to make way for new data). Create 
|	the new block (`newb') before deleting `oldb', then call Reassign to replace `oldb' in `blockList' with `newb'. 
|	Assumes `oldb' exists and is in `blockList'.
*/
void NxsReader::Reassign(
  NxsBlock *oldb,	/* a pointer to the block object soon to be deleted */
  NxsBlock *newb)	/* a pointer to oldb's replacement */
	{
	NxsBlock *prev = NULL;
	NxsBlock *curr = blockList;
	newb->SetNexus(this);

	for (;;)
		{
		if (curr == NULL || curr == oldb)
			break;
		prev = curr;
		curr = curr->next;
		}

	assert(curr != NULL);

	newb->next = curr->next;
	if (prev == NULL) 
		blockList = newb;
	else 
		prev->next = newb;
	curr->next = NULL;
	curr->SetNexus(NULL);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `blockList' data member still equals NULL, returns true; otherwise, returns false. `blockList' will not be equal
|	to NULL if the Add function has been called to add a block object to the list.
*/
bool NxsReader::BlockListEmpty()
	{
	return (blockList == NULL ? true : false);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function was created for purposes of debugging a new NxsBlock. This version does nothing; to create an active
|	DebugReportBlock function, override this version in the derived class and call the Report function of `nexusBlock'.
|	This function is called whenever the main NxsReader Execute function encounters the [&spillall] command comment 
|	between blocks in the data file. The Execute function goes through all blocks and passes them, in turn, to this 
|	DebugReportBlock function so that their contents are displayed. Placing the [&spillall] command comment between
|	different versions of a block allows multiple blocks of the same type to be tested using one long data file. Say 
|	you are interested in testing whether the normal, transpose, and interleave format of a matrix can all be read 
|	correctly. If you put three versions of the block in the data file one after the other, the second one will wipe out
|	the first, and the third one will wipe out the second, unless you have a way to report on each one before the next 
|	one is read. This function provides that ability.
*/
void NxsReader::DebugReportBlock(
  NxsBlock &nexusBlock)	/* the block that should be reported */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(nexusBlock)
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Detaches `oldBlock' from the list of NxsBlock objects growing from `blockList'. If `blockList' itself points to 
|	`oldBlock', this function sets `blockList' to point to `oldBlock->next'. Note: the object pointed to by `oldBlock' 
|	is not deleted, it is simply detached from the linked list. No harm is done in Detaching a block pointer that has 
|	already been detached previously; if `oldBlock' is not found in the block list, Detach simply returns quietly. If 
|	`oldBlock' is found, its SetNexus object is called to set the NxsReader pointer to NULL, indicating that it is no 
|	longer owned by (i.e., attached to) a NxsReader object.
*/
void NxsReader::Detach(
  NxsBlock *oldBlock)	/* a pointer to an existing block object */
	{
	assert(oldBlock != NULL);

	// Return quietly if there are not blocks attached
	//
	if (blockList == NULL)
		return;

	if (blockList == oldBlock) 
		{
		blockList = oldBlock->next;
		oldBlock->SetNexus(NULL);
		}
	else 
		{
		// Bug fix MTH 6/17/2002: old version detached intervening blocks as well
		//
		NxsBlock *curr = blockList;
		for (; curr->next != NULL && curr->next != oldBlock;)
			curr = curr->next;

		// Line below can be uncommented to find cases where Detach function is 
		// called for pointers that are not in the linked list. If line below is
		// uncommented, the part of the descriptive comment that precedes this
		// function about "...simply returns quietly" will be incorrect (at least
		// in the Debugging version of the program where asserts are active).
		//
		//assert(curr->next == oldBlock);

		if (curr->next == oldBlock) 
			{
			curr->next = oldBlock->next;
			oldBlock->SetNexus(NULL);
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reads the NxsReader data file from the input stream provided by `token'. This function is responsible for reading 
|	through the name of a each block. Once it has read a block name, it searches `blockList' for a block object to 
|	handle reading the remainder of the block's contents. The block object is responsible for reading the END or 
|	ENDBLOCK command as well as the trailing semicolon. This function also handles reading comments that are outside 
|	of blocks, as well as the initial "#NEXUS" keyword. The `notifyStartStop' argument is provided in case you do not 
|	wish the ExecuteStart and ExecuteStop functions to be called. These functions are primarily used for creating and 
|	destroying a dialog box to show progress, and nested Execute calls can thus cause problems (e.g., a dialog box is 
|	destroyed when the inner Execute calls ExecuteStop and the outer Execute still expects the dialog box to be 
|	available). Specifying `notifyStartStop' false for all the nested Execute calls thus allows the outermost Execute 
|	call to control creation and destruction of the dialog box.
*/
void NxsReader::Execute(
  NxsToken	&token,				/* the token object used to grab NxsReader tokens */
  bool		notifyStartStop)	/* if true, ExecuteStarting and ExecuteStopping will be called */
	{
	char id_str[256];
	currBlock = NULL;

	bool disabledBlock = false;
	string errormsg;

	try
		{
		token.GetNextToken();
		}
	catch (NxsException & x)
		{
		NexusError(token.GetErrorMessage(), 0, 0, 0, x.cResultCode);
		return;
		}

	if (!token.Equals("#NEXUS"))
		{
		errormsg = "Expecting #NEXUS to be the first token in the file, but found ";
		errormsg << token.GetToken() << " instead";
		NexusError(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn(), kCmdFailedGenerateMessage);
		return;
		}

	if (notifyStartStop)
		ExecuteStarting();

	for (;;)
		{
		token.SetLabileFlagBit(NxsToken::saveCommandComments);
		token.GetNextToken();

		if (token.AtEOF())
			break;

		if (token.Equals("BEGIN"))
			{
			disabledBlock = false;
			token.GetNextToken();

			for (currBlock = blockList; currBlock != NULL; currBlock = currBlock->next)
				{
				if (token.Equals(currBlock->GetID()))
					{
					if (currBlock->IsEnabled()) 
						{
						strcpy(id_str, currBlock->GetID().c_str());
						bool ok_to_read = EnteringBlock(id_str);
						if (!ok_to_read) 
							currBlock = NULL;
						else
							{
							currBlock->Reset();

							// We need to back up currBlock, because the Read statement might trigger
							// a recursive call to Execute (if the block contains instructions to execute 
							// another file, then the same NxsReader object may be used and any member fields (e.g. currBlock)
							//  could be trashed.
							//
							NxsBlock *tempBlock = currBlock;	

							try 
								{
								currBlock->Read(token);
								currBlock = tempBlock;
								}
							catch (NxsException & x) 
								{
								currBlock = tempBlock;
								file_pos pos;
								unsigned line, col;
									// some errors are generated without the offending NxsToken, so they don't put in any file location
									// we'll supply the current location here (should  be better than nothing).
								if (x.pos < 1 && line == 0 && col == 0)
									{
									pos = token.GetFilePosition();
									col = token.GetFileColumn();
									line = token.GetFileLine();
									}
								else
									{
									pos = x.pos;
									col = x.col;
									line = x.line;
									}
								if (currBlock->GetErrorMessage().length() > 0)
									NexusError(currBlock->GetErrorMessage(), pos, line, col, x.cResultCode);
								else
									NexusError(x.msg, pos, line, col, x.cResultCode);
								currBlock = NULL;
								return;
								}	// catch (NxsException x) 
							ExitingBlock(id_str /*currBlock->GetID()*/);
							}	// else
						}	// if (currBlock->IsEnabled()) 

					else
						{
						disabledBlock = true;
						SkippingDisabledBlock(token.GetToken());
						}
					break;
					}	// if (token.Equals(currBlock->GetID()))
				}	// for (currBlock = blockList; currBlock != NULL; currBlock = currBlock->next)

			if (currBlock == NULL)
				{
				string currBlockName = token.GetToken();

				if (!disabledBlock) 
					SkippingBlock(currBlockName);

				for (;;)
					{
					token.GetNextToken();

					if (token.Equals("END") || token.Equals("ENDBLOCK")) 
						{
						token.GetNextToken();

						if (token.GetTokenReference() != ';') 
							{
							errormsg = "Expecting ';' after END or ENDBLOCK command, but found ";
							errormsg << token.GetToken() << " instead";
							NexusError(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn(), kCmdFailedGenerateMessage);
							return;
							}
						break;
						}

					if (token.AtEOF()) 
						{
						errormsg = "Encountered end of file before END or ENDBLOCK in block ";
						errormsg += currBlockName;
						NexusError(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn(), kCmdFailedGenerateMessage);
						return;
						}
					}	// for (;;)
				}	// if (currBlock == NULL)
			currBlock = NULL;
			}	// if (token.Equals("BEGIN"))

		else if (token.Equals("&SHOWALL"))
			{
			for (NxsBlock*  showBlock = blockList; showBlock != NULL; showBlock = showBlock->next)
				{
				DebugReportBlock(*showBlock);
				}
			}

		else if (token.Equals("&LEAVE"))
			{
			break;
			}

		} // for (;;)

	if (notifyStartStop)
		ExecuteStopping();

	currBlock = NULL;
	}
/*----------------------------------------------------------------------------------------------------------------------
|	This function may be used to report progess while reading through a file. For example, the NxsAllelesBlock class 
|	uses this function to report the name of the population it is currently reading so the user doesn't think the 
|	program has hung on large data sets.
*/
void NxsReader::OutputComment(
  const string &comment)	/* a comment to be shown on the output */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(comment)
#	endif
	}
	
#endif

