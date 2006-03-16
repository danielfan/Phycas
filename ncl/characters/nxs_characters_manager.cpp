#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/taxa/nxs_taxa_manager.hpp"	//@POL Mark, moved this line above line that includes nxscharactersmanager.h
#include "ncl/characters/nxs_characters_manager.hpp"
#include "ncl/characters/nxs_characters_block.hpp"
#include <boost/bind.hpp>
#include "ncl/misc/eliminated_index_slider.hpp"
#include "ncl/characters/stored_matrix.hpp"	
#include "ncl/command/nxs_set_cmd_param.hpp"
#include "ncl/output/nxs_output.hpp"
using std::string;

NxsIndexSetCmdOption *InstantiateCharSetOption(const string &n, NxsIndexSet *manip, const string &d, NxsCharactersManager *tmPtr, bool  persist, CmdPermissionLevel pLevel)
	{
	return InstantiateSetOptionGeneric<NxsCharactersManager>(n, manip, d, tmPtr, persist , "Char", pLevel);
	}

/*--------------------------------------------------------------------------------------------------
|	adds another stored matrix to storedMatrices (this matrix will point to the the `matrix' arg
*/
CmdResult NxsCharactersManager::AddCharacters(
  ncl::DiscreteMatrixShPtr matrix, 
  unsigned totalNChar, 
  const NxsIndexSet *eliminated,
  const LabelMap *charLabels, 
  const StateLabelMap *charStateLabels)
  	{
  	unsigned nCharsBefore = 0;
	if (!matrices.empty())
		{
		VecMatrixShPtr::const_iterator mIt = matrices.end();
		--mIt;
		nCharsBefore = (*mIt)->getTotalNChars() + (*mIt)->getNCharsBefore();
		}
	StorageMatrixShPtr newMat;
	newMat = StorageMatrixShPtr(new ncl::StoredMatrix(nCharsBefore, matrix, totalNChar, eliminated, charLabels, charStateLabels));
	matrices.push_back(newMat);
	nChars += totalNChar;
  	AlertListeners(NxsCharListener::kCharsAdded);
  	return kCmdSucceeded;
  	}

CmdResult NxsCharactersManager::NewBlockRead(NxsBlock *newBlock)
	{
	NxsCharactersBlock * newChars = dynamic_cast<NxsCharactersBlock *>(newBlock);
	if (newChars != NULL)
		{
		const unsigned numNewChars = newChars->GetTotalNumCharacters();
		const NxsIndexSet &elim = newChars->GetEliminatedCharacters();
		const ncl::DiscreteMatrixShPtr matrixPtr = newChars->GetMatrix();	//we can alias this because a characters block won't alter the matrix after it is read.
		const LabelMap & charLabels = newChars->GetCharLabels();
		const StateLabelMap &charStateLabels = newChars->GetStateLabels();
		if (this->reportNewCharacters && NxsOutput::GetOutputStreamPtr() != NULL)
			newChars->Report(*NxsOutput::GetOutputStreamPtr());
		CmdResult r =  AddCharacters(matrixPtr, numNewChars, &elim, &charLabels, &charStateLabels);
		newChars->RemoveMatrixAndReset();
		IndicesIncreased(GetSize(), true);
		return r;
		}
	assert(0);
	return kCmdSucceeded;
	}

void NxsCharactersManager::AlertListeners(NxsCharListener::CharChangeType ct)
	{
	std::for_each(listeners.begin(), listeners.end(), boost::bind(&NxsCharListener::CharsChanged, _1, this, ct));
	}

NxsCharactersManager::NxsCharactersManager(NxsTaxaManager & taxaMgr)
  	:NxsTaxaListener(&taxaMgr),
  	nChars(0),
	nxsTaxaMgr(taxaMgr),
#	if defined(CIPRES_MODULE_CONFIG_H)
		reportNewCharacters(false)
#	else
		reportNewCharacters(true)
#	endif
  	{
  	}

NxsCharactersManager::~NxsCharactersManager()
	{
	}
	
NxsCharsBlockAlias NxsCharactersManager::GetCharsBlockReader()
	{
	SharedCharsBlockPtr temp = SharedCharsBlockPtr(new NxsCharactersBlock(nxsTaxaMgr, *this));
	charBlocks.push_back(temp);
	return NxsCharsBlockAlias(temp);	
	}	

	
