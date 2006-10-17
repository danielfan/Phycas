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
#include "phypy/src/ncl/taxa/nxs_alternative_taxa_block.hpp"
#include "phypy/src/ncl/nxs_exception.hpp"
#include "phypy/src/ncl/taxa/nxs_taxa_manager.hpp"	//!oops not part of ncl
using std::string;
/*----------------------------------------------------------------------------------------------------------------------
|	Flushes all locally stored taxa info and gets pointer to the default taxa block from the NxsTaxaManager.
|	Should be called whenever a NxsBlock (derived from NxsAlternativeTaxaBlock) is about to be read.
|	(ie call in the Reset() function or something called by Reset
|	derived class is responsible for resetting dimensionSettings->newTaxa !!! (Note CharactersBlock and DataBlock do this
|	 in CanReadBlockType() )
*/
void NxsAlternativeTaxaBlock::ResetTaxaInfo()
	{
	taxonPos.clear();
	nTaxaInMatrix = 0;
	createdNewTaxa = false;
	advancedCmdRead = false;
	dimensionsRead = false;
	taxLabels.clear();
	labelsRead = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	SHOULD ONLY BE CALLED AFTER StartingCommandThatUsesTaxa when we are creating a new matrix
*/
void NxsAlternativeTaxaBlock::SetTaxonLabel(unsigned ind, const std::string &s)
	{
	NXS_ASSERT(ind < taxLabels.size());
	NXS_ASSERT(IsAddingNewTaxa());
	if (!IsLegalNexusWord(s))
		{
		std::string msg;
		msg << s.c_str() << " is not a legal taxon label."; 
		throw NxsException(s);
		}
	taxLabels.SetLabel(ind, s);
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Called after the DIMENSIONS command is parsed.  
|	Makes sure there are taxa or Newtaxa was set.  If Newtaxa was true a new taxa block is created (via call to 
|	TaxonManager::GetNewTaxaBlockReader()
*/
CmdResult  NxsAlternativeTaxaBlock::HandleDimensions(NxsDimensionsSettings *s)
	{
	dimensionsRead = true;
	if (s->newTaxa)
		{
		if (s->nTaxa == UINT_MAX)
			throw NxsException("NTax was not specified in the DIMENSIONS command");
		createdNewTaxa = true;
 		NXS_ASSERT(taxLabels.empty());
		VecString tmp;
 		FillVectorWithNumbers(tmp, 1, s->nTaxa+1);
		taxLabels.AppendLabels(tmp);
		}
	return FinishHandlingDimensions(s);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new NxsAlternativeTaxaBlock with tb as the source of all taxa info
*/
NxsAlternativeTaxaBlock::NxsAlternativeTaxaBlock(NxsTaxaManager & inTaxaMgr)
  :taxaMgr(inTaxaMgr)
	{
  	ResetTaxaInfo();
  	}
 
void NxsAlternativeTaxaBlock::StartingCommandThatUsesTaxa(const string &cmdName) 
 	{
 	StartingAdvancedCommand();
 	if (!IsAddingNewTaxa() && taxLabels.empty())
 		{
 		dimensionSettings->nTaxa = taxaMgr.GetSize();
 		for (unsigned i = 0; i < dimensionSettings->nTaxa; ++i)
 			taxLabels.push_back(taxaMgr.GetLabel(i));
 		}
 	if (taxLabels.empty())
		{
		string s;
		s << "A " << (IsAddingNewTaxa() ? "NTAX must be read from the DIMENSIONS command" : "TAXA block must be read") << " before the "  << cmdName << " command";
		throw NxsException(s);
		}
	nTaxaInMatrix = GetNumTaxa();	//note nTaxaInMatrix will be wrong if not all taxa appear in the matrix (until SetNewNumTaxaAfterReadingData is called)
 	}
 	
/*----------------------------------------------------------------------------------------------------------------------
|	Should be called when ready to read the data (beginning of MATRIX command parsing) to make the taxonPos vector big 
|	enough.
|	fills the taxonPos vector with default indices (UINT_MAX if the order of taxa isn't known yet)
|	The defaults are changed by calls to SetTaxPos()
*/
void NxsAlternativeTaxaBlock::BuildTaxonPosArray(
  bool matchingTaxaBlock)	/*true if the order must be the same as in the taxa block (e.g. there are no labels on the matrix in a characters block) */
  	{
	taxonPos.clear();
	if (taxLabels.empty())
		StartingAdvancedCommand();
	const unsigned nTaxaRecognized = (unsigned)taxLabels.size();
	for (unsigned i = 0; i < nTaxaRecognized; ++i)
		taxonPos.push_back(matchingTaxaBlock ? i : UINT_MAX);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Alerts the NxsTaxaManager to new taxa (if this block added taxa)
|	should be called when no more changes to taxa names are possible.  As the block exits is fine (from 
|	NxsBlock::EndEncountered() works).
*/
void NxsAlternativeTaxaBlock::FinishedManipulatingTaxa() const
	{
	if (createdNewTaxa && !taxLabels.empty())
		taxaMgr.ReplaceTaxa(taxLabels.GetLabels());
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns current index of a taxon in matrix.  Takes the index of the taxon in the taxa block that was used to read
|	the characters block.  Returns -1 if data was not present for the indicated taxon.
*/
unsigned NxsAlternativeTaxaBlock::GetTaxPos(const string &taxName) const
	{
	unsigned i = FindExternalIndexForTaxon(taxName);
	if (i == UINT_MAX)
		{
		string s;
		s << "The taxon " << taxName << " is unrecognized.";
		throw NxsException(s);
		}
	return GetTaxPos(i);
	}

