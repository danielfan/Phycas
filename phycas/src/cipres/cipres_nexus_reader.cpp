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

// This file is included by read_nexus_i.cpp to avoid a bizarre anonymous namespace multiple definition link error that TL is getting on Mac 10.3.9 (gcc 3.3)
#if defined (INCLUDE_TO_AVOID_LINK_ERROR) || defined(POL_PHYCAS)

//#include "phycas/force_include.h"
#include <fstream>
#include <map>
#include "phycas/src/cipres/cipres_nexus_reader.hpp"
#include "phycas/src/oldphycas/taxa_manager.hpp"
#include "phycas/src/cipres/trees_manager.hpp"
#include "phycas/src/oldphycas/characters_manager.hpp"
#include "phycas/src/ncl/characters/nxs_characters_block.hpp"
#include "phycas/src/ncl/taxa/nxs_taxa_block.hpp"
#include "phycas/src/ncl/trees/nxs_trees_block.hpp"
#include "phycas/src/ncl/misc/string_extensions.hpp"
#include "phycas/src/ncl/nxs_token.hpp"
#include "phycas/src/cipres/CipresDataMatrixHelper.h"
#include "phycas/src/cipres/ConfigDependentHeaders.h"
#include "phycas/src/cipres/AllocateMatrix.hpp"
#include "phycas/src/cipres/CipresAssert.hpp"
#include "phycas/src/cipres/util_copy.hpp"
#include "phycas/src/ncl/characters/stored_matrix.hpp"

#if defined(NO_IDL_TYPES)
#	include "phycas/src/cipres/ReadNexusConstants.h"
#	include "phycas/src/ncl/nxs_exception.hpp"
#else
#	include "CipresIDL/ReadNexusS.h"
	const int NEXUS_TAXA_BLOCK_BIT = CipresIDL::ReadNexus::NEXUS_TAXA_BLOCK_BIT;
	const int NEXUS_TREES_BLOCK_BIT = CipresIDL::ReadNexus::NEXUS_TREES_BLOCK_BIT;
	const int NEXUS_CHARACTERS_BLOCK_BIT = CipresIDL::ReadNexus::NEXUS_CHARACTERS_BLOCK_BIT;
	const int NEXUS_ASSUMPTIONS_BLOCK_BIT = CipresIDL::ReadNexus::NEXUS_ASSUMPTIONS_BLOCK_BIT;
	const int NEXUS_SETS_BLOCK_BIT = CipresIDL::ReadNexus::NEXUS_SETS_BLOCK_BIT;
#endif
using std::map;
using std::ifstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;

const unsigned MAX_CIPR_MATRIX_ELEMENT = (1 << ((sizeof(CIPR_StateSet_t) * 8) - 1)) - 1;

const std::vector<FullTreeDescription> & CipresNexusReader::GetTrees() const
	{
	PHYCAS_ASSERT(phoTreesMgr);
	return phoTreesMgr->GetTrees();
	}

#if defined(POL_PHYCAS) //@POL 10-Nov-2005 added
vector<string> CipresNexusReader::GetTaxLabels() const
	{
	PHYCAS_ASSERT(phoTaxaMgr);
	unsigned ntaxa = phoTaxaMgr->GetNumTaxa();
	vector<string> tax_labels;
	tax_labels.reserve(ntaxa);
	for (unsigned i = 0; i < ntaxa; ++i)
		{
		tax_labels.push_back(phoTaxaMgr->GetUserSuppliedLabel(i));
		}
	return tax_labels;
	}
#endif
	
	// helper function used when we encounter a new state code and the ncl data storage is a single DataStorageType code
CIPR_Datatypes nclToNativeDatatype(const NxsDataType::NxsDatatypesEnum);
typedef map<DataStorageType, CIPR_StateSet_t> NCLToCipresIndexMap;
CIPR_StateSet_t addNewStateSetSingleCode(
  vector<unsigned>	  & tempStateListPos, 
  vector<CIPR_State_t>	  & tempStateList,
  string			  & symbolsForCIPR, 
  const DataStorageType	currNCLState, 
  const NxsDataType	  & datatype,
  NCLToCipresIndexMap &	nclToIndexMap,
  const bool			verbose);

	// helper function used when we encounter a new state code and the ncl data storage is an array of DataStorageType code
typedef vector<const DataStorageType *> VecOfNCLStateArrays;
CIPR_StateSet_t addNewStateSetMultiElementCode(
  vector<unsigned>				  &	tempStateListPos, 
  vector<CIPR_State_t>	  			  &	tempStateList,
  string			  			  &	symbolsForCIPR, 
  const DataStorageType		* const	currNCLState, 
  const NxsDataType	  			  &	datatype,
  VecOfNCLStateArrays 			  &	nclToIndexMap,
  const DataStorageType * const  *&	begVecNCLCodes,
  unsigned  					  &	numNCLCodesStored,
  bool								hasGap,
  const bool						verbose);
  
template<typename T>
unsigned findArrayIndex(unsigned lenToCompare, const T * const * beg, const T * const * const end, const T * const toFind);
void addSymbolForState(string & symbolsForCIPR, const DataStorageType * const  nclState, bool verbose);


string incrementTreesTaxaIndices(const string & zeroBasedNewick)
	{
	string toReturn;
	toReturn.reserve(zeroBasedNewick.length());
	const VecString strVec = NxsToken::TokenizeString(zeroBasedNewick);
	string prevToken;
	for (VecString::const_iterator sIt = strVec.begin(); sIt != strVec.end(); ++sIt)
		{
		unsigned ind;
		if (IsAnUnsigned(*sIt, &ind) && (prevToken == "(" || prevToken == "," || prevToken == ")"))
			toReturn << ind + 1;
		else
			toReturn += *sIt;
		prevToken = *sIt;
		}
	return toReturn;
	}
	
CIPR_Datatypes nclToNativeDatatype(const NxsDataType::NxsDatatypesEnum dt)
	{
	switch (dt)
		{
		case NxsDataType::kDNA:
		case NxsDataType::kNucleotide:
										return CIPR_DNA_Datatype;
		case NxsDataType::kRNA:
										return CIPR_RNA_Datatype;
		case NxsDataType::kProtein:		return CIPR_AA_Datatype;
		case NxsDataType::kStandard:	return CIPR_Generic_Datatype;
		default:
			PHYCAS_ASSERT(0);
		}
		//unreachable
	return CIPR_Generic_Datatype;
	}

void addSymbolForState(string & symbolsForCIPR, const DataStorageType * nclState, const NxsDataType	& datatype, bool verbose)
	{
	try	{
		symbolsForCIPR += datatype.GetStateChar(nclState, false);
#		if !defined(POL_PHYCAS) //POL 29-Oct-2005 to eliminate unwanted output when calling from Python
			if (verbose)
				cout << "symbolsForCIPR.length() = " << (unsigned)symbolsForCIPR.length() <<"\nsymbolsForCIPR = \"" << symbolsForCIPR << "\"\n";
#		endif
		}
	catch (const NxsDataType::XStateNotFound &)
		{
		//@ this need attention, if we see a combination of states without an equate code, 
		// 		what should we push into symbolsForCIPR?
		//	I'm putting in ' ', but PAUP will not be happy if the symbols aren't unique
		symbolsForCIPR += ' ';
		if (verbose)
			cout << "adding space symbolsForCIPR.length() = " << (unsigned)symbolsForCIPR.length() << "\nsymbolsForCIPR = \"" << symbolsForCIPR << "\"\n";
		}
	}
	
CIPR_StateSet_t addNewStateSetSingleCode(
  vector<unsigned>	  & tempStateListPos, 
  vector<CIPR_State_t>	  & tempStateList,
  string			  & symbolsForCIPR, 
  const DataStorageType	currNCLState, 
  const NxsDataType	  & datatype,
  NCLToCipresIndexMap &	nclToIndexMap,
  const bool			verbose)
	{
	const unsigned nStates = datatype.GetNumStates();
	const CIPR_StateSet_t indexToReturn = (CIPR_StateSet_t)tempStateListPos.size();
	tempStateListPos.push_back((unsigned)tempStateList.size());
	vector<CIPR_State_t> newStates;
	for (unsigned i = 0; i < nStates; ++i)
		{
		const DataStorageType nclCode = *datatype.GetStateInBitCode(symbolsForCIPR[i]);
		if (nclCode & currNCLState)
			newStates.push_back(i);
		}
#	if !defined(POL_PHYCAS) //POL 31-Oct-2005 to eliminate unwanted output when calling from Python
		if (verbose)
			cout << "newStates.size() = " << (unsigned)newStates.size() << endl;
#	endif
	tempStateList.push_back((unsigned)newStates.size()); //should be +1, and we should add gaps, but we can't because of the C struct in CipresNativeC.h  
	cip_copy(newStates.begin(), newStates.end(), back_inserter(tempStateList));
	addSymbolForState(symbolsForCIPR, &currNCLState, datatype, verbose);
	nclToIndexMap[currNCLState] = indexToReturn;
	return indexToReturn;
	}

CIPR_StateSet_t addNewStateSetMultiElementCode(
  vector<unsigned>				  &	tempStateListPos, 
  vector<CIPR_State_t>	  			  &	tempStateList,
  string			  			  &	symbolsForCIPR, 
  const DataStorageType		* const	currNCLState, 
  const NxsDataType	  			  &	datatype,
  VecOfNCLStateArrays 			  &	vecNCLCodes,
  const DataStorageType * const  *&	begVecNCLCodes,
  unsigned  					  &	numNCLCodesStored, 
  bool								hasGap,
  const bool						verbose)
	{
	const unsigned nElementsPerCell = datatype.GetNumElementsPerCode();
	if (hasGap)	 // unlike in the single element case, the gap element is not in our lookup.  So we handle gaps here
		{
		unsigned nonZeroEl = 0;
		for (; nonZeroEl != nElementsPerCell; ++nonZeroEl)
			if (currNCLState[nonZeroEl] != 0)
				break;
		if (nonZeroEl == nElementsPerCell)
			return -1;
		}
	const unsigned nStates = datatype.GetNumStates();
	const unsigned indexToReturn = (unsigned)tempStateListPos.size();
	PHYCAS_ASSERT(indexToReturn == vecNCLCodes.size());
	tempStateListPos.push_back((unsigned)tempStateList.size());
	vector<CIPR_State_t> newStates;
	if (hasGap)
		newStates.push_back(-1);
	for (unsigned i = 0; i < nStates; ++i)
		{
		const DataStorageType * nclCode = datatype.GetStateInBitCode(symbolsForCIPR[i]);
		unsigned currEl = 0;
		for (; currEl != nElementsPerCell; ++currEl)
			if ((nclCode[currEl] & currNCLState[currEl]) != nclCode[currEl])
				break;
		if (currEl == nElementsPerCell)
			newStates.push_back(i);
		}

#	if !defined(POL_PHYCAS) //POL 29-Oct-2005 to eliminate unwanted output when calling from Python
		if (verbose)
			cout << "newStates.size() = " << (unsigned)newStates.size() << endl;
#	endif

	tempStateList.push_back((unsigned)newStates.size()); //should be +1, and we should add gaps, but we can't because of the C struct in CipresNativeC.h  
	cip_copy(newStates.begin(), newStates.end(), back_inserter(tempStateList));
	addSymbolForState(symbolsForCIPR, currNCLState, datatype, verbose);
	if ((long) numNCLCodesStored >= (long)MAX_CIPR_MATRIX_ELEMENT)
		{
		const char * msg = "The number of ambiguity codes seen in this data set cannot be held within the current CIPR_Matrix data structure";
#		if defined(NO_IDL_TYPES)
			throw NxsException(msg);
#		else
			throw CipresIDL::NexusException(msg, "", -1, -1);
#		endif
		}
	vecNCLCodes.push_back(currNCLState);
	begVecNCLCodes = &vecNCLCodes[0];
	numNCLCodesStored = (unsigned)vecNCLCodes.size();
	return indexToReturn;
	}
	
template<typename T>
unsigned findArrayIndex(
  unsigned 					lenToCompare, 
  const T * const * 		beg, 
  const T * const * const 	end, 
  const T * const			toFind)
	{
	for (unsigned elementIndex = 0; beg != end; ++elementIndex, ++beg)
		{
		const T * const curr = *beg;
		unsigned i = 0;
		for (; i < lenToCompare; ++i)
			if (curr[i] != toFind[i])
				break;
		if (i == lenToCompare)
			return elementIndex;
		}
	return UINT_MAX;
	}

CIPR_Matrix toCIPRMatrix(const ncl::NxsDiscreteMatrix & inMatrix, const bool verbose)
	{
	CIPR_Matrix toReturn;
	toReturn.stateList = NULL;
	toReturn.stateListPos = NULL;
	toReturn.matrix = NULL;
	toReturn.symbolsList = NULL;
	
	toReturn.nTax = inMatrix.GetNumRows();
	toReturn.nChar = inMatrix.GetNumColumns();
	const NxsDataType &datatype = inMatrix.GetDataType();
	toReturn.datatype = nclToNativeDatatype(datatype.GetDataTypeCategory());
	toReturn.nStates = datatype.GetNumStates();
	const unsigned nStates = toReturn.nStates;
	const char * symCStr = datatype.GetSymbols();
	string symbolsForCIPR(symCStr, symCStr + nStates);
	PHYCAS_ASSERT(symbolsForCIPR.length() == nStates);
	symbolsForCIPR += '?';
	toReturn.matrix = NewTwoDArray<CIPR_StateSet_t>(toReturn.nTax, toReturn.nChar);
	
	vector<CIPR_State_t> tempStateList;
	vector<unsigned> tempStateListPos;

#	if !defined(POL_PHYCAS) //POL 29-Oct-2005 to eliminate unwanted output when calling from Python
		if (verbose)
			{
			cout << "toCIPRMatrix" << std::endl;
			cout << "nTax = "<< toReturn.nTax << '\n';
			cout << "nChar = "<< toReturn.nChar << '\n';
			cout << "nStates = "<< toReturn.nStates << '\n';
			cout << "datatype = "<< toReturn.datatype << '\n';
			}
#	endif
	
	const unsigned nElementsPerCell = datatype.GetNumElementsPerCode();
	if (nElementsPerCell == 1)
		{
		NCLToCipresIndexMap nclToIndex;  // maps ncl code to the index stored in a CIPR_Matrix
		
		CIPR_State_t currStateSetIndex = 0;  // index of the state we are currently describing in  the tempStateList
		
		for (; currStateSetIndex != (CIPR_State_t)nStates; ++currStateSetIndex)
			{
			const DataStorageType nclCode = *datatype.GetStateInBitCode(symbolsForCIPR[currStateSetIndex]);
			nclToIndex[nclCode] = currStateSetIndex;
			tempStateList.push_back(1);
			tempStateList.push_back(currStateSetIndex);
			tempStateListPos.push_back(2*currStateSetIndex);
			}
		const DataStorageType missingNCLCode = *datatype.GetMissingBitCode();
		nclToIndex[missingNCLCode] = currStateSetIndex; // in must check gap status to see if this is ? or N
		tempStateListPos.push_back(2*currStateSetIndex);
		tempStateList.push_back(nStates + 1);
		tempStateList.push_back(-1);
		for (unsigned i = 0; i < nStates; ++i)
			tempStateList.push_back(i);
		
		
		nclToIndex[0] = -1; // in NCL gaps are zeroed cells, gaps do not appear in stateList in a CIPR_Matrix
#		if defined(DONT_DISCRIMINATE_QM_FROM_N)
				//this code does not pay attention to whether gaps are a present in ambiguity codes (effectively they are ignored and all N's go to ?, all paritial ambiguities with gap go to the ungapped states only)
			CIPR_StateSet_t missingNotGapIndex = -1;
			for (unsigned taxNum = 0; taxNum < toReturn.nTax; ++taxNum)
				{
				CIPR_StateSet_t * row = toReturn.matrix[taxNum];
				const DataStorageType * currNCLState = inMatrix.GetState(taxNum, 0);
				NCLToCipresIndexMap::const_iterator currPair;
				const IndSet & gapsInThisRow = inMatrix.GetGapMatrixRow(taxNum);
				const IndSet::const_iterator noGapIt = gapsInThisRow.end();
				if (verbose)
					cout << "Num gaps = "<< gapsInThisRow.size() << endl;
				for (unsigned charNum = 0; charNum < toReturn.nChar; ++charNum, ++currNCLState)
					{
					currPair = nclToIndex.find(*currNCLState);
					if (currPair == nclToIndex.end())
						row[charNum] = addNewStateSetSingleCode(tempStateListPos, tempStateList, symbolsForCIPR, *currNCLState, datatype, nclToIndex,  verbose);
					else
						{
						const CIPR_StateSet_t ind = currPair->second;
						if ((ind != (CIPR_StateSet_t) nStates) || gapsInThisRow.find(charNum) != noGapIt)
							row[charNum] = ind;
						else
							{
							if (missingNotGapIndex > 0)
								row[charNum] = missingNotGapIndex;
							else
								{
								missingNotGapIndex = tempStateList.size();
								row[charNum] = addNewStateSetSingleCode(tempStateListPos, tempStateList, symbolsForCIPR, *currNCLState, datatype, nclToIndex, verbose);
								}
							}
						}
					if (verbose)
						cout << int(row[charNum]) << ' ';
					}
				if (verbose)
					cout << '\n';
				}
#		else
				//this code is semi-incorrect in that it always emits a stateList for (all states+ a gap) AND one for (all states but no gaps).
				// this shouldn't mess any one up (CIPR_Matrix should be tolerant to extra stateList elements)
		
				// Go ahead and add the non-gapped ambiguity code
			++currStateSetIndex;
			tempStateListPos.push_back((unsigned)tempStateList.size()); 
			tempStateList.push_back(nStates);
			for (unsigned i = 0; i < nStates; ++i)
				tempStateList.push_back(i);
			nclToIndex[missingNCLCode] = currStateSetIndex; // replace the code for the ? with N
			addSymbolForState(symbolsForCIPR, &missingNCLCode, datatype, verbose);
				// loop through ignoring gaps, later we'll correct ambiguity codes with gaps (this keeps the gap checking out of this tight loop) 
			for (unsigned taxNum = 0; taxNum < toReturn.nTax; ++taxNum)
				{
				CIPR_StateSet_t * row = toReturn.matrix[taxNum];
				const DataStorageType * currNCLState = inMatrix.GetState(taxNum, 0);
				NCLToCipresIndexMap::const_iterator currPair;
				for (unsigned charNum = 0; charNum < toReturn.nChar; ++charNum, ++currNCLState)
					{
					//std::cout << "taxNum = " << taxNum << ", charNum = " << charNum << ", *currNCLState = " << int(*currNCLState) << std::endl;
					currPair = nclToIndex.find(*currNCLState);
					if (currPair == nclToIndex.end())
						row[charNum] = addNewStateSetSingleCode(tempStateListPos, tempStateList, symbolsForCIPR, *currNCLState, datatype, nclToIndex,  verbose);
					else
						row[charNum] = currPair->second;
					}
				}
				// now we are ready to separate codes like {-AG} from {AG} and ? from N
			std::map<CIPR_StateSet_t, CIPR_StateSet_t> ungappedToGapped;
			ungappedToGapped[nStates + 1] = nStates;  // map N to ?, if the cell has a gap
			ungappedToGapped[CIPR_StateSet_t(-1)] = CIPR_StateSet_t(-1);  // map N to ?, if the cell has a gap
			for (unsigned taxNum = 0; taxNum < toReturn.nTax; ++taxNum)
				{
				CIPR_StateSet_t * row = toReturn.matrix[taxNum];
				const IndSet & gapsInThisRow = inMatrix.GetGapMatrixRow(taxNum);
				const IndSet::const_iterator endGapIt = gapsInThisRow.end();
				for (IndSet::const_iterator gapIt = gapsInThisRow.begin(); gapIt != endGapIt; ++gapIt)
					{
					const unsigned charNum = *gapIt;
					const CIPR_StateSet_t ungappedStateCode = row[charNum];
					//std::cout << "taxNum = " << taxNum << ", gap index = " << charNum << ", ungappedStateCode = " << int(ungappedStateCode) << std::endl;
					const std::map<CIPR_StateSet_t, CIPR_StateSet_t>::const_iterator toGappedIt = ungappedToGapped.find(ungappedStateCode);
					if (toGappedIt == ungappedToGapped.end())
						{
						CIPR_StateSet_t gappedStateCode = (unsigned)tempStateListPos.size(); 
						tempStateListPos.push_back((unsigned)tempStateList.size());
						const unsigned ungappedStateListIndex = tempStateListPos[ungappedStateCode];
						const CIPR_State_t nUngappedStates = tempStateList[ungappedStateListIndex];
						tempStateList.push_back(nUngappedStates + 1);
						tempStateList.push_back(-1);
						for (unsigned ungStInd = 0; ungStInd < (unsigned) nUngappedStates; ++ungStInd)
							tempStateList.push_back(tempStateList[ungappedStateListIndex + ungStInd]);
						ungappedToGapped[ungappedStateCode] = gappedStateCode;
						row[charNum] = gappedStateCode;
						}
					else
						row[charNum] = toGappedIt->second;
					}
				}

#		endif
		}
	else
		{
		PHYCAS_ASSERT(nStates > 0);
		VecOfNCLStateArrays vecNCLCodes;
		CIPR_State_t currStateSetIndex = 0;
		for (; currStateSetIndex != (CIPR_State_t)nStates; ++currStateSetIndex)
			{
			const DataStorageType * const nclCode = datatype.GetStateInBitCode(symbolsForCIPR[currStateSetIndex]);
			vecNCLCodes.push_back(nclCode);
			tempStateList.push_back(1);
			tempStateList.push_back(currStateSetIndex);
			tempStateListPos.push_back(2*currStateSetIndex);
			}
		vecNCLCodes.push_back(datatype.GetMissingBitCode());
		tempStateListPos.push_back(2*currStateSetIndex);
		tempStateList.push_back(nStates + 1);
		tempStateList.push_back(-1);
		for (unsigned i = 0; i < nStates; ++i)
			tempStateList.push_back(i);
		
		CIPR_StateSet_t missingNotGapIndex = -1;
		
		const DataStorageType * const * begVecNCLCodes = &vecNCLCodes[0];
		unsigned numNCLCodesStore = (unsigned)vecNCLCodes.size();
		
		for (unsigned taxNum = 0; taxNum < toReturn.nTax; ++taxNum)
			{
			CIPR_StateSet_t * row = toReturn.matrix[taxNum];
			const DataStorageType * currNCLState = inMatrix.GetState(taxNum, 0);
			const IndSet & gapsInThisRow = inMatrix.GetGapMatrixRow(taxNum);
			const IndSet::const_iterator noGapIt = gapsInThisRow.end();
			if (verbose)
				cout << "Num gaps = "<< (unsigned)gapsInThisRow.size() << endl;
			for (unsigned charNum = 0; charNum < toReturn.nChar; ++charNum, currNCLState += nElementsPerCell)
				{
				unsigned charIndex  = findArrayIndex<DataStorageType>(nElementsPerCell, begVecNCLCodes, begVecNCLCodes + numNCLCodesStore, currNCLState);
				if (charIndex == UINT_MAX)
					charIndex = addNewStateSetMultiElementCode(
							tempStateListPos, tempStateList, symbolsForCIPR, currNCLState, datatype, vecNCLCodes, 
							begVecNCLCodes, numNCLCodesStore, gapsInThisRow.find(charNum) != noGapIt, verbose);
				else
					{
					if (gapsInThisRow.find(charNum) != noGapIt)
						{ //has gap
						if (charIndex < (unsigned) nStates || tempStateList[tempStateListPos[charIndex] + 1] != -1)
							{ // we found a matching code, but it didn't contain a gap.  Search again
							charIndex  = findArrayIndex<DataStorageType>(nElementsPerCell, begVecNCLCodes + 1 + charIndex, begVecNCLCodes + numNCLCodesStore, currNCLState);
							if	(charIndex == UINT_MAX)
								charIndex = addNewStateSetMultiElementCode(tempStateListPos, tempStateList, symbolsForCIPR, currNCLState, datatype, vecNCLCodes, begVecNCLCodes, numNCLCodesStore, true, verbose);
							}
						}
					else
						{
						if (charIndex >= (unsigned) nStates)
							{
							if (charIndex == (unsigned) nStates)
								{
								if (missingNotGapIndex > 0)
									charIndex = missingNotGapIndex;
								else
									{
									missingNotGapIndex = (unsigned)tempStateList.size();
									charIndex = addNewStateSetMultiElementCode(tempStateListPos, tempStateList, symbolsForCIPR, currNCLState, datatype, vecNCLCodes, begVecNCLCodes, numNCLCodesStore, false, verbose);
									}
								}
							else
								{
								if (tempStateList[tempStateListPos[charIndex] + 1] == -1)
									{ // we found a matching code, but it didn't contain a gap.  Search again
									charIndex  = findArrayIndex<DataStorageType>(nElementsPerCell, begVecNCLCodes + 1 + charIndex, begVecNCLCodes + numNCLCodesStore, currNCLState);
									if	(charIndex == UINT_MAX)
										charIndex = addNewStateSetMultiElementCode(tempStateListPos, tempStateList, symbolsForCIPR, currNCLState, datatype, vecNCLCodes, begVecNCLCodes, numNCLCodesStore, false, verbose);
									}
								}
							}
						}
					}
				row[charNum] = charIndex;
				if (verbose)
					cout << int(row[charNum]) << ' ';
				}
			if (verbose)
				cout << '\n';
			}
		}
		
	toReturn.nObservedStateSets = (unsigned)tempStateListPos.size();
	toReturn.stateList = new CIPR_State_t[tempStateList.size()];
	cip_copy(tempStateList.begin(), tempStateList.end(), toReturn.stateList);

#	if !defined(POL_PHYCAS) //POL 29-Oct-2005 to eliminate unwanted output when calling from Python
		if (verbose)
			{
			cout << "toCIPRMatrix toReturn.stateList" << endl;
			for (unsigned i = 0; i < tempStateList.size(); ++i)
				cout << '\t'<< int(toReturn.stateList[i]) << endl;
			}
#	endif

	toReturn.stateListPos = new unsigned[toReturn.nObservedStateSets];
	cip_copy(tempStateListPos.begin(), tempStateListPos.end(), toReturn.stateListPos);
	const unsigned lenSymbols = (unsigned)symbolsForCIPR.length();
	char * bareCharPtr = new char [lenSymbols + 1];
	cip_copy(symbolsForCIPR.begin(), symbolsForCIPR.end(), bareCharPtr);
	bareCharPtr[lenSymbols] = '\0';
	toReturn.symbolsList = bareCharPtr;
	
	return toReturn;
	}

// Warning: need to be sure that you know the CIPR_Matrix was created with toCIPRMatrix to use this
void freeCIPRMatrix(CIPR_Matrix & mat)
	{
	delete [] mat.stateList;
	mat.stateList = NULL;
	delete [] mat.stateListPos;
	mat.stateListPos = NULL;
	DeleteTwoDArray<CIPR_StateSet_t>(mat.matrix);
	delete [] mat.symbolsList;
	mat.symbolsList = NULL;
	}

CipresNative::DiscreteMatrix * createNativeDiscreteMatrix(const CipresNexusReader & nexusReader, unsigned int charBlockIndex)
	{
	boost::shared_ptr<const PhoCharactersManager> charMgr = nexusReader.GetCharactersManager();
	boost::shared_ptr<const ncl::StoredMatrix> rawMatrix = charMgr->GetMatrix((unsigned) charBlockIndex);
	if (rawMatrix && rawMatrix->matrix)
		{
		CIPR_Matrix temp_mat = toCIPRMatrix(*(rawMatrix->matrix.get()), false /*verboseMode*/);
		try
			{
			CipresNative::DiscreteMatrix * cipresNativeMatrix = new CipresNative::DiscreteMatrix(temp_mat);
			freeCIPRMatrix(temp_mat);
			return cipresNativeMatrix;
			}
		catch(...)
			{
			freeCIPRMatrix(temp_mat);
			throw;
			}
		}
	return NULL;
	}

CipresNexusReader::CipresNexusReader(const int blockCode)
	{
	phoTaxaMgr = boost::shared_ptr<PhoTaxaManager>(new PhoTaxaManager());
	phoTaxaMgr->WarnBeforeClearing(false);
	//PhoTreesManager::gTaxaMgr = phoTaxaMgr;
	phoTreesMgr = boost::shared_ptr<PhoTreesManager>(new PhoTreesManager(*phoTaxaMgr.get()));
	phoCharactersMgr = boost::shared_ptr<PhoCharactersManager>(new PhoCharactersManager(*phoTaxaMgr.get()));
	NxsReader::Add(phoTaxaMgr->GetTaxaBlockReader()); // always add taxa block because other blocks rely on taxa
	if (blockCode & (int)NEXUS_TREES_BLOCK_BIT)
		{
		TreesBlockShPtr treesBlock = phoTreesMgr->GetTreesBlockReader();
		treesBlock->AllowImplicitNames();
		NxsReader::Add(treesBlock);
		}
	if (blockCode & (int)NEXUS_CHARACTERS_BLOCK_BIT)
		NxsReader::Add(phoCharactersMgr->GetCharsBlockReader());
	}

void CipresNexusReader::NexusError(const std::string &msg, file_pos pos, unsigned line, unsigned col, CmdResult , NxsBlock* )
	{
	errorMsg = msg;
	filePositionOfError = pos;
	lineOfError = line;
	columnOfError = col;
#	if defined(POL_PHYCAS) && 0
		throw my_exception(msg);
#	endif
	}
	
void CipresNexusReader::ReadNxsFilePath(NxsInFilePath & filePath)
	{
	errorMsg.clear();
	const std::string fn = filePath.GetFullName();
	if (!filePath.Exists())
		{
		errorMsg << "The file \'"<< fn <<"\' does not exist.";
		throw NxsException(errorMsg.c_str());
		}
	if (filePath.IsDirectory())
		{
		errorMsg << "\'"<< fn <<"\' is a directory, not a file.";
		throw NxsException(errorMsg.c_str());
		}
	ifstream inStream;
	if (!filePath.Open(inStream))
		{
		errorMsg << "The file \'"<< fn <<"\' could not be opened.";
		throw NxsException(errorMsg.c_str());
		}
	string msg;
	msg << "Executing " << fn;
	//@	NxsOutputOperationStatusID fileOpID = outputMgrRef.StartStatusDisplay(msg, false);

	fileStack.push(fn);
	NxsToken tokenStream(inStream);
	ReadNxsTokenStream(tokenStream);
	fileStack.pop();
	msg.clear();
	msg << "Finished reading " << fn;
	//@	outputMgrRef.EndStatusDisplay(fileOpID, msg);
	}

void CipresNexusReader::ReadStringAsNexus(const std::string & s)
	{
	NxsToken tokenStream(s);
	ReadNxsTokenStream(tokenStream);
	}	

void CipresNexusReader::ReadNxsTokenStream(NxsToken & tokenStream)
	{
	errorMsg.clear();
	const bool successfullyRead = Execute(tokenStream);
	if (!successfullyRead)
		throw NxsException(errorMsg.c_str(), (long) filePositionOfError, (long) lineOfError, (long) columnOfError);
	}

//INCLUDE_TO_AVOID_LINK_ERROR hack
#endif

