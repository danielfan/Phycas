//#define INCLUDE_TO_AVOID_LINK_ERROR
// This file inclusion avoids a bizarre anonymous namespace multiple definition link error that TL is getting on Mac 10.3.9 (gcc 3.3)
//#warning using macro to include nxs_command_output.cpp TEMPORARY HACK!
//#include "ncl/command/nxs_command_output.cpp"


//#include "phycas/force_include.h"
#include "pyphy/src/ncl/nxs_token.hpp"
#include "pyphy/src/ncl/output/nxs_output.hpp"
#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/characters/nxs_characters_block.hpp"
#include "pyphy/src/ncl/characters/nxs_characters_manager.hpp"
#include "pyphy/src/ncl/taxa/nxs_taxa_manager.hpp"
#include "pyphy/src/ncl/nxs_exception.hpp"
#include "pyphy/src/ncl/nxs_token.hpp"
#include "pyphy/src/ncl/command/nxs_cmd_param.hpp" // needed because we have to check if commands have been read in HandleFormat
#include "pyphy/src/ncl/misc/nxs_data_type_inl.hpp"
#include "pyphy/src/ncl/output/nxs_output_stream.hpp"

using std::map;
using std::vector;
using std::string;
using std::set;
using std::pair;
using ncl::DiscreteMatrixShPtr;
using ncl::NxsDiscreteMatrix;

void  	NxsCharactersBlock::Report(
  NxsOutputStream& outS) const
	{
	outS << '\n' << GetID() << " block contains ";
	if(GetNumTaxaWithData() == 0)
		outS << "no taxa";
	else if(GetNumTaxaWithData() == 1)
		outS << "one taxon";
	else
		outS << GetNumTaxaWithData() << " taxa";
	outS << " and ";
	string s;
	AppendNumberThenWord(s, GetTotalNumCharacters(), "character");
	outS << s << "\n  Data type is \"";
	switch(formatSettings.dataTypeIndex)
		{
		case NxsDataType::kDNA:
			outS << "DNA\"\n";
			break;
		case NxsDataType::kRNA:
			outS << "RNA\"\n";
			break;
		case NxsDataType::kNucleotide:
			outS << "nucleotide\"\n";
			break;
		case NxsDataType::kProtein:
			outS << "protein\"\n";
			break;
		case NxsDataType::kContinuous:
			outS << "continuous\"\n";
			break;
		default:
			outS << "standard\"\n";
		}
	if(formatSettings.respectingCase)
		outS << "  Respecting case\n";
	else
		outS << "  Ignoring case\n";
	if(formatSettings.tokens)
		outS << "  Multicharacter tokens allowed in data matrix\n";
	else
		outS << "  Data matrix entries are expected to be single symbols\n";
	if(formatSettings.labels && formatSettings.transposing)
		outS << "  Character labels are expected on left side of matrix\n";
	else if(formatSettings.labels && !formatSettings.transposing)
		outS << "  Taxon labels are expected on left side of matrix\n";
	else
		outS << "  No labels are expected on left side of matrix\n";
	if(!charLabels.empty())
		{
		outS << "  Character and character state labels:\n";
		for(unsigned k = 0; k < GetTotalNumCharacters(); ++k) 
			{
			const string charKLabel = GetLabel(k);
			outS << '\t' << (1 + k) << '\t' << charKLabel << "\n";

			// output state labels if any are defined for this character
			StateLabelMap::const_iterator cib = charStates.find(k);
			if(cib != charStates.end())
				{
				unsigned ns = (unsigned) (*cib).second.size();
				for(unsigned m = 0; m < ns; ++k)
					outS << "\t\t" << (*cib).second[m] << "\n";
				}
			}
		}
	if(formatSettings.transposing && formatSettings.interleaving)
		outS << "  Matrix transposed and interleaved\n";
	else if(formatSettings.transposing && !formatSettings.interleaving)
		outS << "  Matrix transposed but not interleaved\n";
	else if(!formatSettings.transposing && formatSettings.interleaving)
		outS << "  Matrix interleaved but not transposed\n";
	else
		outS << "  Matrix neither transposed nor interleaved\n";
	outS << "  Missing data symbol is '" << formatSettings.missingSymbol << "\'\n";
	if(formatSettings.matchSymbol != '\0')
		outS << "  Match character is '" << formatSettings.matchSymbol << "\'\n";
	else
		outS << "  No match character specified\n";
	if(formatSettings.gapSymbol != '\0')
		outS << "  Gap character specified is '" << formatSettings.gapSymbol << "\'\n";
	else
		outS << "  No gap character specified\n";
	const char *symbolsList = dataType.GetSymbols();
	outS << "  Valid symbols are: " << symbolsList << "\n";
	const map<char, string>	eqMap = dataType.GetEquates();
	if(!eqMap.empty()) 
		{
		outS << "  Equate macros in effect:\n";
      	for(map<char, string>::const_iterator i = eqMap.begin(); i != eqMap.end(); ++i) 
			outS << "    " << (*i).first << " = " << (*i).second << '\n';
		}
	else
		outS << "  No equate macros have been defined\n";
	if(nEliminated == 0)
		outS << "  No characters were eliminated\n";
	else 
		{
		outS << "  The following characters were eliminated:\n";
		NxsIndexSet::const_iterator k(eliminated.begin());
		for(; k != eliminated.end(); ++k)
			outS << "    " << ((*k)+1) << "\n";
		}
	outS << "  Data matrix:" << "\n";
	DebugShowMatrix(outS, false, "    ");
	
	}
/**
 * @method DebugShowMatrix [int:protected]
 * @param out [ostream&] output stream on which to print matrix
 * @param use_matchchar [bool] if true, matchchar symbol used; otherwise, states shown for all taxa
 * @param marginText [char*] string to print first on each line
 *
 * Provides a dump of the contents of the matrix variable.  Useful for testing
 * whether data is being read as expected.  The default for marginText is NULL,
 * which has the effect of placing the matrix output flush left.  If each line
 * of output should be prefaced with a tab character, specify marginText = "\t".
 */
void NxsCharactersBlock::DebugShowMatrix(
  NxsOutputStream& outS, 
  bool use_matchchar, 
  const char* marginText) const
	{
	const unsigned width = NxsAlternativeTaxaBlock::taxLabels.GetMaxLabelLength();
	unsigned first_taxon = UINT_MAX;
	const unsigned nCharTotal = GetTotalNumCharacters();
	
	for(unsigned j = 0; j < GetNumTaxa(); ++j)
		{
		unsigned taxPos = GetTaxPos(j);
		if (taxPos != UINT_MAX)
			{
			if(marginText != NULL)
				outS << marginText;
			const string &currTaxonLabel = NxsAlternativeTaxaBlock::taxLabels.GetLabel(j);
			outS << currTaxonLabel;
			// print out enough spaces to even up the left edge of the matrix output
			unsigned currTaxonLabelLen = (unsigned) currTaxonLabel.size();
			unsigned diff = width - currTaxonLabelLen;
			for(unsigned k = 0; k < diff + 5; ++k)
				outS << ' ';
			if(first_taxon == UINT_MAX)
				first_taxon = taxPos;
			for(unsigned currChar = 0; currChar < nCharTotal; ++currChar) 
				{
	     		unsigned k = charIndexAdjuster.GetElementIndexFromLocalIndex(currChar);
	     		if(k != UINT_MAX)
		 			ShowStateLabels(outS, taxPos, k, (use_matchchar ? first_taxon : UINT_MAX));
				}
			outS << '\n';
			}
		}
	
	for(unsigned j = 0; j < GetNumTaxa(); ++j)
		{
		if (GetTaxPos(j) == UINT_MAX)
			outS << marginText << "No data for " << NxsAlternativeTaxaBlock::taxLabels.GetLabel(j) << '\n';	
		}
	outS << ncl::endl;
	}


/**
 * @method ShowStateLabels [void:protected]
 * @param out [ostream&] the output stream on which to write
 * @param i [int] the taxon, in range [0..ntax)
 * @param j [int] the character, in range [0..nchar)
 * @param first_taxon [int] the index of the first taxon (if -1, don't use matchchar)
 *
 * Looks up the state(s) at row i, column j of matrix and writes it (or them)
 * to out.  If there is uncertainty or polymorphism, the list of states is
 * surrounded by the appropriate set of symbols (i.e., parentheses for polymorphism,
 * curly brackets for uncertainty).  If 'tokens' is in effect, the output takes
 * the form of the defined state labels; otherwise, the correct symbol is
 * looked up in symbols and output.
 */
void NxsCharactersBlock::ShowStateLabels(
  NxsOutputStream& outS, 
  unsigned rowIndex, unsigned colIndex, unsigned first_taxon) const
	{
	if (rowIndex == UINT_MAX || colIndex == UINT_MAX)
		return;
	if(dataMatrix->IsMissing(rowIndex, colIndex))
		{
		outS << formatSettings.missingSymbol;
		return;
		}
	else if(dataMatrix->IsGap(rowIndex, colIndex))
		{
		outS << formatSettings.gapSymbol;
		return;
		}
	const DataStorageType *stateCode = dataMatrix->GetState(rowIndex, colIndex);
	if(formatSettings.tokens)
		{
		unsigned n = dataType.CountStates(stateCode);
		if(n == 1) 
			{
			bool use_matchchar = false;
			if(first_taxon != UINT_MAX && rowIndex != first_taxon) 
				{
				const DataStorageType * firsts = dataMatrix->GetState(first_taxon, colIndex);
				if(dataType.StatesComp(firsts, stateCode) == 0)
					use_matchchar = true;
				}
			if(use_matchchar)
				outS << formatSettings.matchSymbol;
			else 
				{
				StateLabelMap::const_iterator ci = charStates.find(colIndex);
				// OPEN ISSUE: need to eliminate state labels for characters that have
				// been eliminated
				bool printed = false;
				if (ci != charStates.end())
					{
					vector<unsigned> ords = dataType.GetOrdinationsOfStates(stateCode);
					if (ci->second[ords[0]] != " ")
						{
						printed = true;
						outS << "  " << ci->second[ords[0]];
						}
					}
				if (!printed)
					outS << "  " << dataType.GetStateChar(stateCode) << "[<-no label found]";
				}
			}
		else 
			{
			vector<unsigned> ords = dataType.GetOrdinationsOfStates(stateCode);
			
					//TODO: handle matchchar possibility here too
			//
			if(dataMatrix->IsPolymorphic(rowIndex, colIndex))
				outS << "  (";
			else
				outS << "  {";
			for(vector<unsigned>::iterator s = ords.begin(); s != ords.end(); ++s) 
				{
				StateLabelMap::const_iterator ci = charStates.find(colIndex);
				if(ci == charStates.end() || (*ci).second[*s] == " ")
					outS << "  " << dataType.GetStateCharFromIndex(*s) << "[<-no label found]";
				else		// show label at index number s in LabelList at ci
					outS << "  " << (*ci).second[*s];
				}
			if(dataMatrix->IsPolymorphic(rowIndex, colIndex))
				outS << ')';
			else
				outS << '}';
			}
		}
	else 
		{
		if(first_taxon != UINT_MAX && rowIndex > first_taxon && (dataType.StatesComp(stateCode, dataMatrix->GetState(first_taxon, colIndex)) == 0))
			outS << '.';
		else 
			outS << dataType.GetStatesAsNexusString(stateCode, dataMatrix->IsPolymorphic(rowIndex,colIndex), false);
		}
	}


CmdResult NxsCharactersBlock::EndEncountered()
	{
	NxsAlternativeTaxaBlock::FinishedManipulatingTaxa();
	return charactersMgr.NewBlockRead(this);
	}

void NxsCharactersBlock::ResetDataType(NxsDataType::NxsDatatypesEnum dataTypeIndex)
	{
	dataType = NxsDataType(dataTypeIndex, dataTypeIndex == NxsDataType::kStandard);
	//std::cout <<" end of ResetDataType" << std::endl; //@POL 27-Oct-2005 looked like debugging code, so I commented it out
	}


CmdResult  NxsCharactersBlock::HandleFormat(NxsFormatCmdSettings *s)
	{
	formatSettings = *s;
	if(formatSettings.dataTypeIndex == NxsDataType::kContinuous) 
		{
		errorMsg = "Sorry, continuous character matrices have not yet been implemented";
		throw NxsException(errorMsg);
		}
	
	//	check all of the difficult inter-dependencies in the format command
	//
	if (!IsLegalNexusChar(formatSettings.missingSymbol))
		ThrowIllegalNexusChar(formatSettings.missingSymbol, "MISSING");
	if (formatSettings.gapSymbol != '\0')
		{
		if (!IsLegalNexusChar(formatSettings.gapSymbol))
			ThrowIllegalNexusChar(formatSettings.gapSymbol, "GAP");
		if (formatSettings.gapSymbol == formatSettings.missingSymbol)
			throw NxsException("the GAP character and MISSING character cannot be identical");
		if (formatSettings.gapSymbol == formatSettings.matchSymbol)
			throw NxsException("the GAP character and MATCHCHAR cannot be identical");
		}
	if (formatSettings.matchSymbol != '\0')
		{
		if (!IsLegalNexusChar(formatSettings.matchSymbol))
			ThrowIllegalNexusChar(formatSettings.matchSymbol, "MATCHCHAR");
		if (formatSettings.matchSymbol == formatSettings.missingSymbol)
			throw NxsException("the missing character and MATCHCHAR character cannot be identical");
		}
	bool allowRC = false;
	if (formatSettings.dataTypeIndex == NxsDataType::kContinuous)
		{
		if (formatSettings.rawSymbols.length() > 0)
			throw NxsException("the SYMBOLS keyword is not allowed if the DATATYPE is CONTINUOUS");
		if (!statesFormatCI.get()->HasBeenRead())
			{
			//	Individuals is the default for the states format command if datatype = continuous
			//
			formatSettings.stateFormIndex = 1; // depends on "STATESPRESENT|INDIVIDUALS|COUNT|FREQUENCY" order of choices
			formatSettings.stateFormName = "INDIVIDUALS";
			}
		if (!formatSettings.tokens)
			{
			if (tokensCI.get()->HasBeenRead())
				throw NxsException("TOKENS must be used when the DATATYPE is CONTINUOUS");
			formatSettings.tokens = false;
			}
		}
	else if (formatSettings.dataTypeIndex == NxsDataType::kDNA || formatSettings.dataTypeIndex == NxsDataType::kRNA || formatSettings.dataTypeIndex == NxsDataType::kNucleotide)
		{
		if (formatSettings.tokens)
			throw NxsException("TOKENS cannot be used with NUCLEOTIDE data");
		}
	else
		allowRC = true;
	if (!allowRC && formatSettings.respectingCase)
		throw NxsException("the RESPECTCASE option can only be used with the STANDARD DATATYPE");
	//	get the symbols list ready for reading in the matrix
	//
	ResetDataType((NxsDataType::NxsDatatypesEnum) formatSettings.dataTypeIndex);
	dataType.SetGapChar('\0');
	dataType.SetMatchChar('\0');
	dataType.SetMissingChar('\0');
	
	if (formatSettings.rawSymbols.length() > 0)
		{
		if (formatSettings.dataTypeIndex == NxsDataType::kStandard)
			dataType.ReplaceSymbols(formatSettings.rawSymbols, formatSettings.respectingCase);
		else
			dataType.AddSymbols(formatSettings.rawSymbols, formatSettings.respectingCase);
		}
	//	missing, match and gap characters can't be in the symbols list.
	//
	const string fullSymbolsList = string(dataType.GetSymbols());
	string errCharName;
	char errorChar  = '\0';
	if (fullSymbolsList.find(formatSettings.gapSymbol) != string::npos)
		{
		errCharName = "GAP";
		errorChar = formatSettings.gapSymbol;
		}
	else if (fullSymbolsList.find(formatSettings.matchSymbol) != string::npos)
		{
		errCharName = "MATCH";
		errorChar = formatSettings.matchSymbol;
		}
	else if (fullSymbolsList.find(formatSettings.missingSymbol) != string::npos)
		{
		errCharName = "MISSING";
		errorChar = formatSettings.missingSymbol;
		}
	if (!errCharName.empty())
		{
		errorMsg << "The " << errCharName << " character(" << errorChar << ") cannot be identical to a state's symbol";
		throw NxsException(errorMsg);
		}
		
	dataType.SetGapChar(formatSettings.gapSymbol);
	dataType.SetMatchChar(formatSettings.matchSymbol);
	dataType.SetMissingChar(formatSettings.missingSymbol);
			
	if (formatSettings.rawEquate.length() > 0)
		dataType.AddEquates(formatSettings.rawEquate);
	return kCmdSucceeded;
	}


/**
 * @method HandleTokenState [int:protected]
 * @param token [NxsToken&] the token used to read from in
 * @param j [int] the character index, in range [0..nchar)
 * @throws NxsException
 *
 * Called from HandleNextState to read in the next state when 'tokens' is in effect.
 * Looks up state in character states listed for the character to make
 * sure it is a valid state, and returns state's value (0, 1, 2, ...).
 * Note: does NOT handle adding the state's value to matrix. Save the return
 * value, let's call it k, and use the following command to add it to matrix:
 * matrix->AddState(i, j, k);
 */
unsigned NxsCharactersBlock::ReadTokenState(
  NxsToken	&token, 
  unsigned origInd) const
	{
	// token should be one of the character states listed for character origInd
	// in charStates
	//
   	StateLabelMap::const_iterator bagIter = charStates.find(origInd);
	if(bagIter == charStates.end()) 
   		{
		errorMsg << "No states were defined for character " << (1 + origInd);
		throw NxsException(errorMsg, token);
		}
	const VecString &statesVec(bagIter->second);
   // TO DO: this section is very UGLY - need to find some cleaner way of comparing
   // the token string to the strings representing valid characters states
   // in the LabelList associated with character j
   //
	VecString::const_iterator cit;
	if(formatSettings.respectingCase)
		cit = find(statesVec.begin(), statesVec.end(), token.GetTokenReference());
	else
		{
		NStrCaseInsensitiveEquals compF(token.GetTokenReference());
		cit = find_if(statesVec.begin(), statesVec.end(), compF);
		}
		
	if(cit == statesVec.end()) 
		{
		errorMsg << "Character state " << token.GetTokenReference() << " not defined for character " << (1 + origInd);
		throw NxsException(errorMsg, token);
		}
	// ok, the state has been identified, so return the state's internal
	// representation. That is, if the list of state labels was
	// "small medium large" and "small" was specified in the data file,
	// state saved in matrix would be 0 (it would be 1 if "medium" were
	// specified in the data file, and 2 if "large" were specified in the
	// data file).
	return (unsigned) distance(statesVec.begin(), cit);
	}



bool NxsCharactersBlock::ReadNextStateToken(
  NxsToken& token, 
  unsigned i, 
  unsigned j, 
  unsigned origInd)
	{
	NXS_ASSERT(formatSettings.tokens);
	++token;
	// handle case in which TOKENS was specified in the FORMAT command
	// token should be in one of the following forms: "{"  "a"  "bb"
	//
	bool polymorphism = (token.GetTokenReference() == '(');
	bool uncertainty  = (token.GetTokenReference() == '{');

	if(!uncertainty && !polymorphism) 
		{
		unsigned k = ReadTokenState(token, origInd);
		dataMatrix->SetStateIndex(i, j, k);
		}
	else 
		{
     	bool prevTokenWasTilde = false;
     	unsigned first = UINT_MAX;
     	unsigned last;
     	dataMatrix->ZeroState(i,j);
		for(;;)
			{
			// OPEN ISSUE: What about newlines if interleaving? I'm assuming
			// that the newline must come between characters to count.
			++token;
			if (polymorphism && token.GetTokenReference() == ')') 
				{
				if(prevTokenWasTilde) 
					{
					errorMsg = "Range of states still being specified when ')' encountered";
					throw NxsException(errorMsg, token);
					}
				break;
				}
			else if (uncertainty && token.GetTokenReference() == '}') 
				{
				if(prevTokenWasTilde) 
					{
					errorMsg = "Range of states still being specified when '}' encountered";
					throw NxsException(errorMsg, token);
					}
				break;
				}
			else if (token.GetTokenReference() == '~') 
				{
				if(first == UINT_MAX) 
					{
					errorMsg = "Tilde character ('~') cannot precede token indicating beginning of range";
					throw NxsException(errorMsg, token);
					}
				prevTokenWasTilde = true;
				}
			else if(prevTokenWasTilde) 
				{
				// Add all states from first+1 to last, then reset prevTokenWasTilde to 0
           		//
				last = ReadTokenState(token, origInd);
				if(last <= first) 
					{
					errorMsg << "Last state in specified range (" << token.GetTokenReference() << ") must be greater than the first";
					throw NxsException(errorMsg, token);
					}

				for(unsigned k = first+1; k <= last; ++k)
           			dataMatrix->AddStateIndex(i, j, k);
				prevTokenWasTilde = false;
           		first = UINT_MAX;
				}
			else 
				{
	           // Add current state, then set first to that state's value
	           // State's value is its position within the list of states
	           // for that character
	           //
				first = ReadTokenState(token, origInd);
				dataMatrix->AddStateIndex(i, j, first);
				}
			}
		if(polymorphism)
			dataMatrix->SetPolymorphic(i, j);
		}
	return true;
	}


bool NxsCharactersBlock::ReadTransposedMatrix(NxsToken& token)
	{
	const unsigned nCharTotal = GetTotalNumCharacters();
	const unsigned totalNTax = GetNumTaxa();
	const bool readLabels = formatSettings.labels;
	const bool interleaving = formatSettings.interleaving;
	const bool tokens = formatSettings.tokens;
	unsigned currTaxon = UINT_MAX;
	unsigned firstTaxonInPage = 0;
	unsigned lastTaxonInPage = 	totalNTax;
	bool newCharLabels = (needCharLabels && charLabels.empty());
	vector<unsigned> indToCharPosMap;
	charIndexAdjuster.BuildOldToNewIndexMap(&indToCharPosMap, nCharTotal);
	
	
		// if currTaxon equals nTaxInMatrix, then we've just finished reading the last
		// interleave page and thus should break from the outer loop
		// Note that if we are not interleaving, this will still work since
		// lastTaxon is initialized to ntaxTotal and never changed
		//
	for (unsigned page = 0; currTaxon != totalNTax; ++page)
		{
		for(unsigned currChar = 0; currChar < nCharTotal; ++currChar)
			{
			if(readLabels)
				{
				// this should be the character label
				//
				++token;
				unsigned charInd = FindIndexFromCharLabels(token.GetTokenReference());
				bool isANumber;
				if (page == 0 && newCharLabels)
					{
					if (charInd != UINT_MAX) 
						{
						errorMsg << "Data for this character (" << token.GetTokenReference() << ") has already been saved";
						throw NxsException(errorMsg, token);
						}

					// Labels provided, need to add them to charLabels list.
					// saving character labels even for characters that have been eliminated.
					//
					if (!IsALegalNexusLabelForObjectN(token.GetTokenReference(), currChar, &isANumber))
						{
						errorMsg << token.GetTokenReference() << " is an illegal " << GetCharLabelsCmdName() << " for " << GetDatumName() << " number " << (currChar+1) << '.';
						if (isANumber)
							errorMsg << "  If a number is used as a" << GetCharLabelsCmdName() << ", it must identical to the number of that " << GetDatumName();
						throw NxsException(errorMsg, token);
						}
					if (!isANumber)
						charLabels[currChar] = token.GetTokenReference();
					}
				else 
					{
					// either not first interleaved page or character labels previously defined
					//
					if(charInd == UINT_MAX) 
						{
						if (!IsALegalNexusLabelForObjectN(token.GetTokenReference(), currChar, &isANumber))
							{
							errorMsg << token.GetTokenReference() << " is an illegal " << GetCharLabelsCmdName() << " for " << GetDatumName() << " number " << (currChar+1) << '.';
							throw NxsException(errorMsg, token);
							}
						if (!isANumber)
							{
							errorMsg << "Could not find character named " << token.GetTokenReference() << " among stored character labels";
							throw NxsException(errorMsg, token);
							}
						charInd = currChar;
						}
					else if (charInd != currChar)
						{
						// make sure user has not duplicated data for a single character or
						// changed the order in which characters appear in different interleave
						// pages
						//
						if(page == 0)
							errorMsg << "Data for this character (" << token.GetTokenReference() << ") has already been saved";
						else
							errorMsg = "Ordering of characters must be identical to that in first interleave page";
						throw NxsException(errorMsg, token);
						}
					}
				} 
			
			//************************************************
			//******** Beginning of loop through taxa ********
			//************************************************
			if (interleaving)
				token.AlterTokenReading(NxsToken::kNewlineIsToken);
			try	
				{
				// 	 it is possible that character currChar has been ELIMINATEd, we need to keep track 
				//	 of  the positon in characters matrix by iterating posInStored
				//	
				unsigned charPosInMatrix = indToCharPosMap[currChar];
				NxsDiscreteMatrix *bareDataMatrix = dataMatrix.get();
				for(currTaxon = firstTaxonInPage; currTaxon < lastTaxonInPage; ++currTaxon)
					{
					// ok will be 0 only if a newline character is encountered before
					// taxon i processed
					//
					bool ok ;
					if (!tokens)
						{
						token.ReadSingleCharacter();
						NXS_ASSERT(token.GetTokenLength() > 0);
						if (token.GetTokenReference()[0] == '\n')
							ok = false;
						else
							{
							ok = true;
							dataType.ReadNextState<NxsDiscreteMatrix>(bareDataMatrix, token, currTaxon, charPosInMatrix);
							}
						}
					else
						ok = ReadNextStateToken(token, currTaxon, charPosInMatrix, currChar);
							
					if(!ok)
						{
						NXS_ASSERT(interleaving);
						if (lastTaxonInPage < totalNTax && currTaxon != lastTaxonInPage) 
							{
							errorMsg = "Each line within an interleave page must comprise the same number of taxa";
							throw NxsException(errorMsg, token);
							}

						lastTaxonInPage = currTaxon;
						}
					} // innermost loop (over taxa)
				}
			catch (...)
				{
				token.AlterTokenReading(NxsToken::kNewlineIsNotToken);
				throw;
				}
			if (interleaving)
				token.AlterTokenReading(NxsToken::kNewlineIsNotToken);
			} // middle loop (over characters)
		firstTaxonInPage = lastTaxonInPage;
		lastTaxonInPage = totalNTax;
		} // outer loop (over interleave pages)
	++token;
	token.ThrowIfNot(";");
	return true;
	}


bool NxsCharactersBlock::ReadStdMatrix(
  NxsToken& token)
	{
	const unsigned nCharTotal = GetTotalNumCharacters();
	unsigned maxNTax = NxsAlternativeTaxaBlock::GetNumTaxa();
	const bool readLabels = formatSettings.labels;
	const bool interleaving = formatSettings.interleaving;
	const bool tokens = formatSettings.tokens;
	bool labelIsNumber;
	unsigned currChar = 0;
	unsigned firstCharInCurrPage = 0;
	bool needNextTaxonName = true; 	// used in case we accidentally read a taxon name from the next interleaved page (will only happen when interleaving, without newtaxa and when not all of the taxa from the taxa block are present)
	vector<unsigned> indToCharPosMap;
	charIndexAdjuster.BuildOldToNewIndexMap(&indToCharPosMap, nCharTotal);
	vector<unsigned>::const_iterator thisPagesFirstPosInMatrix = indToCharPosMap.begin();
	vector<unsigned>::const_iterator currPosInMatrix = indToCharPosMap.begin();
		// if currChar equals ncharTotal, then we've just finished reading the last
		// interleave page and thus should break from the outer loop
		// Note that if we are not interleaving, this will still work since
		// 		lastCharInCurrPage is initialized to ncharTotal and never changed 
		//		so currChar == nCharTotal when it exits the loop over taxa
		//
	for(unsigned page = 0; currChar != nCharTotal; ++page)
		{
		unsigned lastCharInCurrPage = nCharTotal;
		for(unsigned rowIndex = 0; rowIndex < maxNTax; ++rowIndex)
			{
			unsigned positionInTaxaBlock = rowIndex;
			if(readLabels)
				{
				if (needNextTaxonName)
					++token;	// This should be the taxon label
				
				if (token.GetTokenReference() == ';' && lastCharInCurrPage == nCharTotal && !IsAddingNewTaxa() && page == 0)
					{
					//	we hit a semicolon before readin all of the taxa from the taxa block.  This is tolerated as long as the user didn't specify ntax
					//	and we aren't partially through reading the data for the taxa that are included.
					//
					SetNewNumTaxaAfterReadingData(rowIndex);
					return true;
					}
				needNextTaxonName = true;
				if (page == 0 && IsAddingNewTaxa())
					{
					//	the label supplied should be unique (or a number identical to the taxon number)
					//
					unsigned prevInd = GetTaxPosDontThrow(token.GetTokenReference());
					if (prevInd != UINT_MAX && prevInd != rowIndex)
						{
						errorMsg << "Data for this taxon (" << token.GetTokenReference() << ") has already been saved";
						throw NxsException(errorMsg, token);
						}
					// replacing the default name (#) with a name
					//
					if (!IsALegalNexusLabelForObjectN(token.GetTokenReference(), rowIndex, &labelIsNumber))
						{
						errorMsg << token.GetTokenReference() << " is not a legal taxon name.";
						if (labelIsNumber)
							errorMsg << "  If a number is used as a taxon label, it must identical to the number of that taxon in the matrix.";
						throw NxsException(errorMsg, token);
						}
					SetTaxonLabel(rowIndex, token.GetTokenReference());	
					}
				else
					{
					// Cannot assume taxon in same position in
					// taxa block. Set up taxonPos array so that we can look up
					// the correct row in matrix for any given taxon
					//
					positionInTaxaBlock = FindExternalIndexForTaxon(token.GetTokenReference());
					if (positionInTaxaBlock == UINT_MAX) 
						{
						if (token.GetTokenReference() == ';' && rowIndex == 0)
							errorMsg << "Unexpected ; (after only " << currChar << " characters were read)";
						else
							errorMsg << "Could not find taxon named " << token.GetTokenReference() << " among stored taxon labels";
						throw NxsException(errorMsg, token);
						}
					if (page == 0)
						{
						if (GetTaxPos(positionInTaxaBlock) != UINT_MAX)
							{
							if (interleaving && GetTaxPos(positionInTaxaBlock) == 0)
								{
								// we've repeated a taxon name without reaching the rowIndex == maxNTax.  
								// 	This is only allowed if we are interleaving (and he user didn't specify ntax)
								//
								needNextTaxonName = false;	//flags the fact that we've already read the taxon name
								maxNTax = rowIndex;
								SetNewNumTaxaAfterReadingData(rowIndex);
								break;	// breaks to the loop over interleave pages
								}
							else
								{
								errorMsg << "Data for this taxon (" << token.GetTokenReference() << ") has already been saved";
								throw NxsException(errorMsg, token);
								}
							}
						}
					else
						{
						if(GetTaxPos(positionInTaxaBlock) != rowIndex) 
							{
							errorMsg = "Ordering of taxa must be identical to that in first interleave page";
							throw NxsException(errorMsg, token);
							}
						}
					}
				if (page == 0)
					SetTaxPosInMatrix(positionInTaxaBlock, rowIndex);
				}
						
			//******************************************************
			//******** Beginning of loop through characters ********
			//******************************************************
			if (interleaving)
				token.AlterTokenReading(NxsToken::kNewlineIsToken);
			try	
				{
				NxsDiscreteMatrix *bareDataMatrix = dataMatrix.get();
				currPosInMatrix = thisPagesFirstPosInMatrix;
				for(currChar = firstCharInCurrPage; currChar < lastCharInCurrPage; ++currChar, ++currPosInMatrix)
					{
					// 	Because some characters might be eliminated, posInStored does not 
					//	necessarily equal currChar
					//
					bool ok;
					
					if (!tokens)
						{
						token.ReadSingleCharacter();
						NXS_ASSERT(token.GetTokenLength() > 0);
						if (token.GetTokenReference()[0] != '\n')
							{
							ok = true;
							dataType.ReadNextState<NxsDiscreteMatrix>(bareDataMatrix, token, rowIndex, *currPosInMatrix);
							}
						else
							ok = false;
						}
					else
						ok = ReadNextStateToken(token, rowIndex, *currPosInMatrix, currChar);
						
					if(!ok)
	            		{
	            		// ok will be false only if a newline character is encountered
						//
						NXS_ASSERT(interleaving);
	            		if(lastCharInCurrPage < nCharTotal && currChar != lastCharInCurrPage) 
							{
							// lastCharInCurrPage ==  nCharTotal for the first taxon in a page.
							//
							errorMsg = "Each line within an interleave page must comprise the same number of characters";
							throw NxsException(errorMsg, token);
							}
						lastCharInCurrPage = currChar;
						}
					} // innermost loop (over characters)
				}
			catch (...)
				{
				token.AlterTokenReading(NxsToken::kNewlineIsNotToken);
				throw;
				}
			if (interleaving)
				token.AlterTokenReading(NxsToken::kNewlineIsNotToken);
			} // middle loop (over taxa)
		//	to read the next page, advance charPosIt and firstCharInCurrPage
		//
		thisPagesFirstPosInMatrix =  currPosInMatrix;
		firstCharInCurrPage = lastCharInCurrPage;
		} // outer loop (over interleave pages)
	++token;
	token.ThrowIfNot(";");
	return true;
	}


bool  NxsCharactersBlock::ParseMatrix(
  NxsToken & token)
	{
	StartingCommandThatUsesTaxa("Matrix");
	unsigned expectedNTax = GetNumTaxa();
	if (expectedNTax == 0) 
		{
		errorMsg << "Must precede " << GetID() << " block with a TAXA block or specify NEWTAXA and NTAX in the DIMENSIONS command";
		throw NxsException(errorMsg, token);
		}
	
	dataMatrix = DiscreteMatrixShPtr(new NxsDiscreteMatrix(expectedNTax, nChar, dataType));
	//	if we're reading a non-transposed matrix with labels then we need to build the taxonPos array while reading the matrix
	//	if not the ordering has to be identical to the taxablock ordering
	//
	BuildTaxonPosArray(!(formatSettings.labels && !(formatSettings.transposing)));
	
	bool matRead;
	if(formatSettings.transposing)
		matRead = ReadTransposedMatrix(token);
	else
		matRead = ReadStdMatrix(token);
	NXS_ASSERT(matRead);
	//  now the ordering of the taxa in the matrix is known, so we alert the AlternativeTaxaBlock interface
	//	that we are done introducing new taxa, or deciding on the order of taxa 
	//
	FinishedManipulatingTaxonPos();
	if (GetNumTaxaWithData() < expectedNTax)
		{
		for (unsigned i = expectedNTax -1; i >= GetNumTaxaWithData(); --i)
			dataMatrix->RemoveRow(i);
		dataMatrix->TrimExcessCapacity();
		}
	return kCmdSucceeded;
	}

CmdResult  NxsCharactersBlock::HandleEliminate(NxsEliminateCmdSettings *s)
	{
	eliminated = s->toEliminate;
	nEliminated = eliminated.size();
	charIndexAdjuster.SetEliminated(eliminated, GetTotalNumCharacters());
	nChar -= nEliminated;
	return kCmdSucceeded;
	}


void NxsCharactersBlock::ReadCharLabel(
  NxsToken &token,
  set<string> &uniqNames,
  unsigned index,
  bool save)
  	{
	if(token.GetTokenReference() != ' ')
		{
		string s(token.GetTokenReference());
		if (!formatSettings.respectingCase)
			Capitalize(s);
		pair< set<string>::iterator ,bool> ret = uniqNames.insert(s);
		if (!ret.second)
			{
			errorMsg << GetCharLabelsCmdName() << "s must be unique (" << s << " was repeated)";
			throw NxsException(errorMsg, token);
			}
		}
	bool isANumber;
	if (!IsALegalNexusLabelForObjectN(token.GetTokenReference(), index, &isANumber))
		{
		errorMsg << token.GetTokenReference() << " is an illegal " << GetCharLabelsCmdName() << " for " << GetDatumName() << " number " << (index+1) << '.';
		if (isANumber)
			errorMsg << "  If a number is used as a" << GetCharLabelsCmdName() << ", it must identical to the number of that " << GetDatumName();
		throw NxsException(errorMsg, token);
		}
	if(save && !isANumber )
		charLabels[index] = token.GetTokenReference();
	++token;
	}

bool  NxsCharactersBlock::ParseCharLabels( 
  NxsToken &token)
	{
	if (NoAdvancedCommandsHaveBeenRead())
		StartingAdvancedCommand();
	if (!charLabels.empty())
		{
		errorMsg << "The " << GetCharLabelsCmdName() << "s command cannot follow a " << GetCharStateLabelsCmdName() << "s command.";
		throw NxsException(errorMsg, token);
		}
		
	set<string> uniqNames;
	++token;
	for(unsigned index = 0; token.GetTokenReference() != ';'; ++index)
		{
		if (index >= GetTotalNumCharacters())
			{
			errorMsg << "The number of " << GetCharLabelsCmdName() << "s supplied cannot exceed " << GetTotalNumCharacters() << " (the " << GetNumCharsName() << " specified in the DIMENSIONS command)";
			throw NxsException(errorMsg, token);
			}
		bool save = !IsEliminated(index);
		ReadCharLabel(token, uniqNames, index, save);
		} 
	return true;
	}

//	return indicates whether or not a ; was encountered
//
bool NxsCharactersBlock::ReadStateLabels(
  NxsToken &token, 
  unsigned index, 
  unsigned nStatesInSymbols, 
  bool save)
	{
	LabelList stList;
	set<string> stSet;
	bool semiColon = false;
	const bool respectCase = formatSettings.respectingCase;
	for(unsigned x = 1; token.GetTokenReference() != ','; ++x)
		{
		if (x > nStatesInSymbols)
			{
			errorMsg << "There only " << nStatesInSymbols << " states according to the FORMAT command, but at least " << x << " state labels were supplied for " << GetDatumName() << ' ' << (index+1) << " in the " << GetCharStateLabelsCmdName() << " command.";
			throw NxsException(errorMsg, token);
			}
		stList.push_back(token.GetTokenReference());
		if (token.GetTokenReference() != ' ')
			{
			string s(token.GetTokenReference());
			if (!respectCase)
				Capitalize(s);
			pair< set<string>::iterator ,bool> ret = stSet.insert(s);
			if (!ret.second)
				{
				errorMsg << "The " << GetStateLabelsCmdName() << " for a " << GetDatumName() << " must be unique (" << s << " was repeated).";
				throw NxsException(errorMsg, token);
				}
			}
		++token;
		if (token.GetTokenReference() == ';') 
			{
			semiColon = true;
			break;
			}
  		} 
	if (save)
		charStates[index] = stList;
	++token;
	return semiColon;
	}


bool  NxsCharactersBlock::ParseStateLabels( 
  NxsToken &token)
	{
	if (NoAdvancedCommandsHaveBeenRead())
		StartingAdvancedCommand();
	if (!charStates.empty())
		{
		errorMsg << "The " << GetStateLabelsCmdName() << " command cannot follow a " << GetCharStateLabelsCmdName() << " command.";
		throw NxsException(errorMsg, token);
		}
	const unsigned nStatesInSymbols = dataType.GetNumStates();
	++token;
	for(; token.GetTokenReference() != ';';)
		{
		const unsigned n = ReadCharacterIndex(token);
		if (ReadStateLabels(token, n, nStatesInSymbols, !IsEliminated(n)))
			break;
	  	} 
	needCharLabels = false;
	return true;
	}
	
unsigned NxsCharactersBlock::ReadCharacterIndex(
  NxsToken &token)	
	{
		// token should be the character number; create a new association
	UInt u;
	if (!IsAnUnsigned(token.GetTokenReference(), &u))
		{
		errorMsg << "Expecting a " << GetDatumName() << " number, but found " << token.GetTokenReference() << " in the " << GetCharStateLabelsCmdName() << "command.";
		throw NxsException(errorMsg, token);
		}
	if (u >  GetTotalNumCharacters() || u < 1) 
		{
		errorMsg << "Invalid " << GetDatumName() << " number (" << token.GetTokenReference() << ") found in " << GetCharStateLabelsCmdName() << " command (";
		if (u <  GetTotalNumCharacters())
			errorMsg << "greater than the " << GetNumCharsName() << " specified in the DIMENSIONS command)";
		else
			errorMsg << GetDatumName() << " numbers must be positive)";
		throw NxsException(errorMsg, token);
		}
	
	++token;
	return (u-1);
	}
	
bool  NxsCharactersBlock::ParseCharStateLabels( 
  NxsToken &token)
	{
	if (NoAdvancedCommandsHaveBeenRead())
		StartingAdvancedCommand();

	if (!charStates.empty() || !charLabels.empty())
		{
		errorMsg << "The " << GetCharStateLabelsCmdName() << " command cannot follow a " << GetCharLabelsCmdName() << " or a " << GetStateLabelsCmdName() << " command.";
		throw NxsException(errorMsg, token);
		}
		
	set<string> uniqNames;
	const unsigned nStatesInSymbols = dataType.GetNumStates();
	++token;
	for(;;)
		{
		if (token.GetTokenReference() == ';')
			break;
		unsigned n = ReadCharacterIndex(token);
		bool save = !IsEliminated(n);
		if (token.GetTokenReference() != '/') 
			ReadCharLabel(token, uniqNames, n, save);
		if (token.GetTokenReference() == '/') 
			{
			++token;
			if (ReadStateLabels(token, n, nStatesInSymbols, save))
				{
				needCharLabels = false;
		  		return true;
		  		}
		  	}
		else if (token.GetTokenReference() == ';')
	  		{
	  		needCharLabels = false;
	  		return true;
	  		}
		else if (token.GetTokenReference() == ',') 
			++token;
		else
			{
			errorMsg << "Expecting a comma or semicolon here, but found (" << token.GetTokenReference() << ") instead";
			throw NxsException(errorMsg, token);
			}
		}
	return true; 
	}



void NxsCharactersBlock::ResetCmdMgrNxsBlock()
	{
	dataType = NxsDataType(NxsDataType::kStandard, false);
	charLabels.clear();	/* vector of the known taxon labels */
	charStates.clear();
	dataMatrix = DiscreteMatrixShPtr((NxsDiscreteMatrix *)NULL);
	nChar = 0;
	needCharLabels = true;
	NxsAlternativeTaxaBlock::ResetTaxaInfo();
	eliminated.clear();
	nEliminated = 0;
	charIndexAdjuster.clear();
	}
	


/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of characters in the longest taxon name (or number)
*/
unsigned NxsCharactersBlock::GetMaxLabelLength() const
	{
	string s;
	s <<  GetTotalNumCharacters();
	unsigned maxN = (unsigned)s.length();
	for (LabelMap::const_iterator labIt = charLabels.begin(); labIt != charLabels.end(); ++labIt)
		{
		if (labIt->second.length() > maxN)
			maxN = (unsigned)labIt->second.length();
		}
	return maxN;
	}


CmdResult NxsCharactersBlock::FinishHandlingDimensions(NxsDimensionsSettings *s)
	{
	nChar = s->secondDimension;
	NXS_ASSERT(charLabels.empty() && charStates.empty());
	return kCmdSucceeded;
	}
	
bool NxsCharactersBlock::CanReadBlockType(
  const string &s)
	{
	if (EqualsCaseInsensitive(s, "CHARACTERS"))
		{
		SetID("CHARACTERS");
		NxsAlternativeTaxaBlock::dimensionSettings->newTaxa = false;
		NxsAlternativeTaxaBlock::dimensionSettings->nTaxa = NxsAlternativeTaxaBlock::taxaMgr.GetSize();
		}
	else if (EqualsCaseInsensitive(s, "DATA"))
		{
		SetID("DATA");
		NxsAlternativeTaxaBlock::dimensionSettings->newTaxa = true;
		}
	else 
		return false;
	return true;
	}



string NxsCharactersBlock::GetAdvancedCommandName() const
	{
	string s;
	s << "a MATRIX or " << GetCharStateLabelsCmdName();
	return s;
	}


NxsCharactersBlock::NxsCharactersBlock(NxsTaxaManager & inTaxaMgr, NxsCharactersManager & inCharMgr)
	:NxsCommandManagerBlock("CHARACTERS"), 
	NxsAlternativeTaxaBlock(inTaxaMgr),
	dataMatrix(),
	dataType(NxsDataType::kStandard, false),
	charactersMgr(inCharMgr)
	{
	ConstructorInitialization();	// call this virtual function which is overridden in AllelesBlock to avoid prematurely calling InitializeRecognizedCommands before AllelesBlock is constructed
	}

void NxsCharactersBlock::ConstructorInitialization()
	{
	InitializeRecognizedCommands();
	Reset();
	}
