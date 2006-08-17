#ifndef NCL_NXSALTERNATIVETAXABLOCK_H
#define NCL_NXSALTERNATIVETAXABLOCK_H

#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/command/nxs_command_decl.hpp"
#include "pyphy/src/ncl/nxs_basic_manager.hpp"
#include "pyphy/src/ncl/taxa/nxs_taxa_label_reader.hpp"

class NxsDimensionsSettings
	{
	public:
		bool newTaxa;
		unsigned nTaxa;
		unsigned secondDimension;	//e.g. nchar
		unsigned thirdDimension; // alleles block uses three dimensions.
	};
typedef boost::shared_ptr<NxsDimensionsSettings> DimensionsSettingsShPtr;
class NxsTaxaManager;
/*----------------------------------------------------------------------------------------------------------------------
|	Several Blocks depend on a list of taxa, and can be used to introduce new taxa labels (e.g. Characters, data, 
|	alleles, distances, and unaligned blocks).  
|	All of these blocks are derived from NxsAlternativeTaxaBlock so that they use the same interactions with the NxsTaxaManager
|	(and don't repeat code).
|	The class:
|		Creates DIMENSIONS, TAXLABELS commands.  
|		Interacts with the NxsTaxaManager to get or create new taxa block as required
|		Handles taxa block index to local index mapping (GetTaxPos())
|
|	The following assumptions are made about the order of functions being called after construction of the object:
|		Reset() is called when the block is entered (this can be done by calling NxsAlternativeTaxaBlock::Reset from the 
|			derived class Reset() -assuming the derived class is a NxsBlock class).  This lets sets the NxsTaxaBlock
|			pointer to the default taxa block.
|		DIMENSIONS and TAXLABELS commands must precede any "Advanced" Commands (those that need to know the 
|			taxlabels - e.g. MATRIX)
|		There should be one command that creates the TaxPos ordering for the block (typically MATRIX) within this
|			command:  	StartingAdvancedCommand() and BuildTaxonPosArray() are called before the data are read
|						ChangeTaxonLabel() and SetTaxPos() are used to set local taxon pos and change names of taxa (if newtaxa labels are added their initial names are simply their numbers)
|						SetNewNumTaxaAfterReadingData() must be called if all of the taxa in the taxa block do not appear in the data
|						FinishedManipulatingTaxonPos() is called after the taxon positions are set (e.g. end of matrix command)
|		FinishedManipulatingTaxa() is called as the block is left.  This lets the NxsAlternativeTaxaBlock notify the taxamanager if new taxa have been introduced
*/
class NxsAlternativeTaxaBlock: public NxsTaxaLabelReader
	{
	public:
		NxsAlternativeTaxaBlock(NxsTaxaManager & taxaMgr);
		virtual ~NxsAlternativeTaxaBlock(){}
		
		bool				DimensionsHasBeenRead() const;
		unsigned 			GetNumTaxaWithData() const;
		unsigned			GetTaxPos(unsigned i) const;
		unsigned			GetTaxPos(const std::string &taxName) const;
		unsigned			GetTaxPosDontThrow(const std::string &taxName) const;
		
		bool 				NoAdvancedCommandsHaveBeenRead() 	{return !advancedCmdRead;}
		
	protected:
				void				StartingAdvancedCommand();
				void				StartingCommandThatUsesTaxa(const std::string &);
				bool				CreatedNewTaxa() {return createdNewTaxa;}
				unsigned 			GetNumTaxa() const;
		virtual std::string			GetAdvancedCommandName() const	{ return "MATRIX";}
		virtual std::string			GetNewTaxaName() const {return "NewTaxa";}
		virtual std::string			GetNumTaxaName() const {return "NTax";}
 		virtual bool 	  			GetNewTaxaDefault() const {return false;}
 
				void				BuildTaxonPosArray(bool matchingTaxaBlock);
				void				SetTaxonLabel(unsigned ind, const std::string &n);
				NxsCommandShPtr	CreateUnfinishedLinkCommand();
				NxsCommandShPtr	CreateUnfinishedDimensionsCommand(NxsTestShPtr);
				NxsCommandShPtr	GetTaxLabelsCommand();
				
				unsigned 			FindExternalIndexForTaxon(const std::string & s) const;

		virtual CmdResult			FinishHandlingDimensions(NxsDimensionsSettings *) = 0;
				void				FinishedManipulatingTaxa() const;
				void				FinishedManipulatingTaxonPos()	{}//currently unused hook to be called when taxon positions in the local block are known.  May never be needed.
				bool				IsAddingNewTaxa();	//should be const, but to do so we need a const test interface
				CmdResult			HandleDimensions(NxsDimensionsSettings *);
				void				ResetTaxaInfo();
				void				SetNewNumTaxaAfterReadingData(unsigned i);
				void				SetTaxPosInMatrix(unsigned taxBlockInd, unsigned localInd);
			
		DimensionsSettingsShPtr dimensionSettings;
	private:

		std::vector<unsigned> 	taxonPos;	 /* UINT_MAX means eliminated, maps NxsTaxaBlock's index to index in matrix */
		unsigned				nTaxaInMatrix;
		bool					createdNewTaxa;
		bool					advancedCmdRead;
		bool					dimensionsRead;
	protected:
		NxsTaxaManager		  & taxaMgr;
	};
	
/*----------------------------------------------------------------------------------------------------------------------
|	Sets advancedCmdRead and gets all of the known taxa from the taxaMgr (if new taxa weren't read)
*/
inline bool  NxsAlternativeTaxaBlock::DimensionsHasBeenRead() const
	{
	return dimensionsRead;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Sets advancedCmdRead and gets all of the known taxa from the taxaMgr (if new taxa weren't read)
*/
inline void NxsAlternativeTaxaBlock::StartingAdvancedCommand() 
 	{
 	advancedCmdRead = true;
 	}
 	
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the block is reading new taxa (as opposed to relying on the taxamanager)
*/
inline bool NxsAlternativeTaxaBlock::IsAddingNewTaxa() 
	{
	return (dimensionSettings->newTaxa && !NxsTaxaLabelReader::labelsRead);
	}	

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current number of taxa active in the block.  Before SetNewNumTaxaAfterReadingData() is called (before
|	MATRIX has been read) this will be 0.
*/
inline unsigned NxsAlternativeTaxaBlock::GetNumTaxaWithData() const
	{
	return nTaxaInMatrix;
	}
	
inline void NxsAlternativeTaxaBlock::SetNewNumTaxaAfterReadingData(
  unsigned i)
	{
	NXS_ASSERT( i <= nTaxaInMatrix);
	nTaxaInMatrix = i;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current index of a taxon in matrix.  Takes the index of the taxon in the taxa block that was used to read
|	the characters block.  Returns UINT_MAX if data was not present for the indicated taxon.
*/
inline unsigned NxsAlternativeTaxaBlock::GetTaxPos(unsigned i) const
	{
	if (i < taxonPos.size())
		return taxonPos[i];
	return UINT_MAX;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	specifies that the taxa block's taxaBlockPos-th taxon is in position localPos
*/
inline void NxsAlternativeTaxaBlock::SetTaxPosInMatrix(
  unsigned taxaBlockPos,
  unsigned localPos)
	{
	NXS_ASSERT(taxaBlockPos < taxonPos.size());
	NXS_ASSERT(taxonPos[taxaBlockPos] == UINT_MAX);
	taxonPos[taxaBlockPos] = localPos;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns current index of a taxon in matrix.  Takes the index of the taxon in the taxa block that was used to read
|	the characters block.  Returns -1 if data was not present for the indicated taxon.
*/
inline unsigned NxsAlternativeTaxaBlock::GetTaxPosDontThrow(const std::string &taxName) const
	{
	unsigned i = FindExternalIndexForTaxon(taxName);
	return i != UINT_MAX ? GetTaxPos(i) : UINT_MAX;
	}


inline unsigned NxsAlternativeTaxaBlock::FindExternalIndexForTaxon(const std::string &taxName) const
	{
	return taxLabels.FindIndex(taxName, false, false);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	SHOULD ONLY BE CALLED AFTER StartingCommandThatUsesTaxa
*/
inline unsigned NxsAlternativeTaxaBlock::GetNumTaxa() const
	{
	return (unsigned)taxLabels.size();
	}

#endif

