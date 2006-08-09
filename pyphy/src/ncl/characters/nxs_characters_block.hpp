#ifndef NCL_NXSCHARACTERSBLOCK_H
#define NCL_NXSCHARACTERSBLOCK_H

#include "ncl/nxs_defs.hpp"
#include "ncl/taxa/nxs_alternative_taxa_block.hpp"
#include "ncl/nxs_basic_manager.hpp"
#include "ncl/nxs_block.hpp"
#include "ncl/misc/nxs_discrete_matrix.hpp"
#include "ncl/misc/nxs_data_type.hpp"
#include "ncl/characters/nxs_characters_block_cmd_settings.hpp"
#include "ncl/misc/eliminated_index_slider.hpp"

class NxsCharactersManager;
class NxsFormatCmdSettings;
class NxsEliminateCmdSettings;
typedef boost::shared_ptr<NxsFormatCmdSettings> CharFormatSettingsPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Class that reads taxa blocks and stores the list of taxon labels.
*/
class NxsCharactersBlock : public NxsCommandManagerBlock , public NxsIndexInfo, public NxsAlternativeTaxaBlock
	{
	public:
				NxsCharactersBlock(NxsTaxaManager & taxaMgr, NxsCharactersManager & charMgr);
	
				//	Accessors
				//
				bool				CanReadBlockType(const std::string &s);
				unsigned			FindIndex(const std::string &label ) const;
				unsigned			FindIndexFromNamesOnly(const std::string &label ) const;
				const NxsIndexSet & GetEliminatedCharacters() const;
				unsigned			GetMaxLabelLength() const;
				unsigned			GetNumCharacters() const;
				unsigned			GetTotalNumCharacters() const;
				unsigned			GetSize() const;
				const std::string & GetUserSuppliedLabel( unsigned i ) const;
				bool				IsInterleave() const;
				bool				IsLabels() const;
				bool				IsRespectCase() const;
				bool				IsTokens() const;
				bool				IsTranspose() const;
				
				void				Report( NxsOutputStream& out ) const;
				//	Modifiers
				//
				void				Clear();
				void				ResetCmdMgrNxsBlock();
				
				//	don't call GetSet	 only added for parsing the eliminate command.
				// 
				const NxsIndexSet *GetSet(const std::string &) const {return NULL;}
				
	protected:
				void			AddRecognizedCommands();
				bool			AddDimensionsCommand(NxsTestShPtr);	
				bool 			AddEliminateCommand(NxsTestShPtr noAdvancedCmdsTest, NxsTestShPtr somCharsExpectedTest);
		virtual bool 			AddFormatCommand(NxsTestShPtr);
 				CmdResult		EndEncountered();
				bool 			ParseTaxLabels(NxsToken &);	
				bool 			ParseDimensions(NxsToken &);	
				
	private:
				void			DebugShowMatrix(NxsOutputStream& out, bool, const char*) const;
				void 			ShowStateLabels(NxsOutputStream& outS, unsigned rowIndex, unsigned colIndex, unsigned first_taxon) const;
				void			BuildCharPosArray();
				CmdResult		HandleEliminate(NxsEliminateCmdSettings *);
				CmdResult		HandleFormat(NxsFormatCmdSettings *);
				CmdResult		FinishHandlingDimensions(NxsDimensionsSettings *);
				bool 			ParseCharLabels(NxsToken & );
				bool 			ParseCharStateLabels(NxsToken & );
				bool 			ParseStateLabels(NxsToken & );
		virtual bool  			ParseMatrix(NxsToken & );
		
				bool			IsEliminated(unsigned ind) const;
				unsigned		FindIndexFromCharLabels(const std::string &label) const;
				
				void 			ReadCharLabel(NxsToken &token, std::set<std::string> &uniqNames, unsigned index, bool save);
				unsigned 		ReadCharacterIndex(NxsToken &token);
				bool	 		ReadNextStateToken(NxsToken& token, unsigned i, unsigned j, unsigned origInd);
				bool 			ReadStateLabels(NxsToken &token, unsigned index, unsigned nStatesInSymbols,  bool save);
				bool  			ReadStdMatrix(NxsToken & );
				unsigned		ReadTokenState(NxsToken	&token, unsigned origInd) const;
				bool  			ReadTransposedMatrix(NxsToken & );
				void 			ResetDataType(NxsDataType::NxsDatatypesEnum);
 				void			SetIndexOffset(unsigned i);
 		//	Functions overridden in AllelesBlock so that CharLabels, CharStateLabels and StateLabels code can be shared
 		virtual void 			ConstructorInitialization();
 		virtual std::string		GetCharLabelsCmdName() const {return "CharLabels";}
 		virtual std::string		GetCharStateLabelsCmdName() const {return "CharStateLabels";}
 		virtual std::string		GetStateLabelsCmdName() const {return "StateLabels";}
		virtual std::string		GetAdvancedCommandName() const;

 		virtual	std::string		GetNumCharsName() const {return "NChar";}
		virtual	std::string		GetDatumName(bool plural = false) const { return plural ? "characters" : "character";}
 		virtual bool 	  		GetNewTaxaDefault() const;

 		enum NxsNucleotideStateCode
				{
				kNucA = 0, 
				kNucC = 1, 
				kNucG = 2,
				kNucT = 3
				};
		enum NxsAminoAcidAtateCode
				{
				kAA_A = 0, kAA_C = 1, kAA_D = 2, kAA_E = 3, kAA_F = 4, kAA_G = 5, kAA_H = 6, kAA_I = 7, kAA_K = 8, kAA_L = 9, 
				kAA_M = 10, kAA_N = 11, kAA_P = 12, kAA_Q = 13, kAA_R = 14, kAA_S = 15, kAA_T = 16, kAA_V = 17, kAA_W = 18, kAA_Y = 19, 
				kAA_Stop = 20
				};
		enum NxsItemTypes 
				{
				kItem_States = 0, 
				kItem_Min, 
				kItem_Max, 
				kItem_Median, 
				kItem_Average, 
				kItem_Variance, 
				kItem_StdError, 
				kItem_SampleSize
				};
		enum NxsStateFormatsChoices 
				{
				kStFormat_StatesPresent = 0, 
				kStFormat_Individ, 
				kStFormat_Count, 
				kStFormat_Freq
				};

		
		NxsCmdOptionShPtr	statesFormatCI;
		NxsCmdOptionShPtr	tokensCI;
		
		LabelMap 				charLabels;	/* vector of the known taxon labels */
		StateLabelMap			charStates;
		ncl::DiscreteMatrixShPtr		dataMatrix;
		unsigned 				nChar;
		bool 					needCharLabels;
		NxsFormatCmdSettings	formatSettings;
		mutable std::string		scratchSpace;
		NxsDataType				dataType;
		NxsIndexSet				eliminated;
		unsigned				nEliminated;
		EliminatedIndexSlider	charIndexAdjuster;
		NxsCharactersManager  & charactersMgr;
		ncl::DiscreteMatrixShPtr			RemoveMatrixAndReset();
		const ncl::DiscreteMatrixShPtr	GetMatrix() const {return dataMatrix;}
		const LabelMap 			  & GetCharLabels() const {return charLabels;}
		const StateLabelMap		  & GetStateLabels() const {return charStates;}
		friend class NxsCharactersManager;
	};
	
#if defined (NCL_SUPPORT_OLD_NAMES)
	typedef NxsCharactersBlock CharactersBlock;
#endif
/*----------------------------------------------------------------------------------------------------------------------
|	Returns index set of elinated characters GetTotalNumCharacters() gives the max number of characters in this set.
*/
inline const NxsIndexSet &NxsCharactersBlock::GetEliminatedCharacters() const
	{
	return eliminated;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a shared pointer to the dataMatrix and resets the block. 
*/
inline ncl::DiscreteMatrixShPtr NxsCharactersBlock::RemoveMatrixAndReset() 
	{
	ncl::DiscreteMatrixShPtr temp = dataMatrix;
	Reset();
	return temp;
	}

				
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of non eliminated characters;
*/
inline unsigned NxsCharactersBlock::GetNumCharacters() const
	{
	return nChar;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number ofcharacters (including eliminated characters)
*/
inline unsigned NxsCharactersBlock::GetTotalNumCharacters() const
	{
	return nChar + nEliminated;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	returns the number of taxa
*/
inline unsigned NxsCharactersBlock::GetSize() const
	{
	return GetNumCharacters();
	}

				
inline void NxsCharactersBlock::Clear()
	{
	Reset();
	}
	
inline bool NxsCharactersBlock::GetNewTaxaDefault() const
	{
	return (blockID == "DATA");
	}
	

inline bool NxsCharactersBlock::IsEliminated(unsigned ind) const
	{
	return eliminated.IsAMember(ind);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Only checks the charLabels (returns UINT_MAX if label is not stored
*/
inline unsigned NxsCharactersBlock::FindIndexFromCharLabels(
  const std::string &label) const
	{
	NStrCaseInsensitiveEquals compF(label);
	for (LabelMap::const_iterator labIt = charLabels.begin(); labIt != charLabels.end(); ++labIt)
		{
		if (compF(labIt->second))
			return labIt->first;
		}
	return UINT_MAX;
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the index for the taxon label supplied (returns UINT_MAX if the label is not recognized).
*/
inline unsigned NxsCharactersBlock::FindIndex(
  const std::string &label) const
	{
	unsigned n = FindIndexFromCharLabels(label);
	if (n != UINT_MAX)
		return n;
	return GetIndexByTreatingLabelAsNumber(label, GetTotalNumCharacters());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the index for the taxon label supplied (returns UINT_MAX if the label is not recognized).
*/
inline unsigned NxsCharactersBlock::FindIndexFromNamesOnly(
  const std::string &label) const
	{
	return FindIndexFromCharLabels(label);
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Returns the label of taxon i (or throws an exception
*/
inline const std::string & NxsCharactersBlock::GetUserSuppliedLabel(
  unsigned i) const
	{
	return GetLabelFromMap(charLabels, i);
	}



#endif
