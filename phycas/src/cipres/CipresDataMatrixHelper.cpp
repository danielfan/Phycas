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

#include "phycas/src/cipres/CipresDataMatrixHelper.h"
#include "phycas/src/cipres/util_copy.hpp"

using CipresNative::DiscreteMatrix;
using std::string;
using std::vector;
using std::cout;
using std::endl;


	// helper function used when we encounter a new state code and the ncl data storage is a single DataStorageType code
CIPR_Datatypes nclToNativeDatatype(const NxsCharactersBlock::DataTypesEnum);

CIPR_Datatypes nclToNativeDatatype(const NxsCharactersBlock::DataTypesEnum dt)
	{
	switch (dt)
		{
		case NxsCharactersBlock::dna:
		case NxsCharactersBlock::nucleotide:
										return CIPR_DNA_Datatype;
		case NxsCharactersBlock::rna:
										return CIPR_RNA_Datatype;
		case NxsCharactersBlock::protein:		return CIPR_AA_Datatype;
		case NxsCharactersBlock::standard:	return CIPR_Generic_Datatype;
		default:
			PHYCAS_ASSERT(0);
		}
		//unreachable
	return CIPR_Generic_Datatype;
	}

DiscreteMatrix::DiscreteMatrix(const NxsCharactersBlock & cb, bool gapsToMissing)
{

	std::vector<const NxsDiscreteDatatypeMapper *> mappers = cb.GetAllDatatypeMappers();
	if (mappers.size() > 1)
		throw NxsException("too many mappers");
	if (mappers.empty() || mappers[0] == NULL)
		throw NxsException("no mappers");
		
	const NxsDiscreteDatatypeMapper & mapper = *mappers[0];
	const NxsDiscreteStateMatrix & rawMatrix = cb.GetRawDiscreteMatrixRef();

	this->nativeCMatrix.datatype = nclToNativeDatatype(mapper.GetDatatype());
	this->nativeCMatrix.nStates = mapper.GetNumStates();
	const std::string fundamentalSymbols = mapper.GetSymbols();
	const std::string fundamentalSymbolsPlusGaps = mapper.GetSymbolsWithGapChar();
	const bool hadGaps = !(fundamentalSymbols == fundamentalSymbolsPlusGaps);
	
	this->symbolsStringAlias = fundamentalSymbols;
	char missingSym = cb.GetMissingSymbol();	
	CIPR_State_t newMissingStateCode = this->nativeCMatrix.nStates;
	PHYCAS_ASSERT((int)NXS_MISSING_CODE < 0);
	PHYCAS_ASSERT((int)NXS_GAP_STATE_CODE < 0);
	const int sclOffset = std::min( (int)NXS_GAP_STATE_CODE, (int)NXS_MISSING_CODE);
	const int negSCLOffset = -sclOffset;
	const unsigned nMapperStateCodes = mapper.GetNumStateCodes();
	const unsigned recodeVecLen = nMapperStateCodes;
	const unsigned nMapperPosStateCodes = nMapperStateCodes - (hadGaps ? 2 : 1);
	std::vector<CIPR_State_t> recodeVec(recodeVecLen, -2);
	CIPR_State_t * recodeArr = &recodeVec[negSCLOffset];

	if (fundamentalSymbols.length() < this->nativeCMatrix.nStates)
		throw NxsException("Fundamental states missing from the symbols string");
	for (CIPR_State_t i = 0; i < (CIPR_State_t)this->nativeCMatrix.nStates; ++i)
		{
		PHYCAS_ASSERT(mapper.PositionInSymbols(fundamentalSymbols[i]) == (int) i);		
#		if defined (NDEBUG)
			const std::set<int>	 & ss =  mapper.GetStateSetForCode(i);
			PHYCAS_ASSERT(ss.size() == 1);
			PHYCAS_ASSERT(*ss.begin() == i);
#		endif
		stateListAlias.push_back(1);
		stateListAlias.push_back(i);
		stateListPosAlias.push_back((unsigned) 2*i);
		recodeArr[i] = i;
		}

	//NXS_INVALID_STATE_CODE

	recodeArr[NXS_GAP_STATE_CODE] = ((hadGaps && gapsToMissing)? newMissingStateCode : -1);
		
	if (missingSym == '\0')
		missingSym = (hadGaps ? mapper.GetGapSymbol() : '?');
	else
		{
		PHYCAS_ASSERT(NXS_MISSING_CODE == mapper.GetStateCodeStored(missingSym));
		}
	recodeArr[NXS_MISSING_CODE] = newMissingStateCode;
	this->symbolsStringAlias.append(1, missingSym);
	const unsigned nCodesInMissing  = this->nativeCMatrix.nStates + (gapsToMissing ?  0 : 1);
	stateListPosAlias.push_back(2*this->nativeCMatrix.nStates);
	stateListAlias.push_back(nCodesInMissing);
	if (!gapsToMissing)
		stateListAlias.push_back(-1);
	for (CIPR_State_t i = 0; i < (CIPR_State_t)this->nativeCMatrix.nStates; ++i)
		stateListAlias.push_back(i);
	
	CIPR_State_t nextStateCode = newMissingStateCode + 1;
	for (int i = (int)this->nativeCMatrix.nStates; i < (int) nMapperPosStateCodes; ++i)
		{
		const std::set<int>	 &ss = mapper.GetStateSetForCode( i);
		const unsigned ns = ss.size();
		const bool mapToMissing  = (!mapper.IsPolymorphic(i) && (nCodesInMissing + 1 == ns || nCodesInMissing == ns));
		if (mapToMissing)
			recodeArr[i] = newMissingStateCode;
		else
			{
			recodeArr[i] = nextStateCode++;
			stateListPosAlias.push_back(stateListAlias.size());
			stateListAlias.push_back(ns);
			for (std::set<int>::const_iterator sIt = ss.begin(); sIt != ss.end(); ++sIt)
				stateListAlias.push_back((CIPR_State_t) *sIt);
			std::string stateName = mapper.StateCodeToNexusString(i);
			if (stateName.length() != 1)
				this->symbolsStringAlias.append(1, ' ');
			else
				this->symbolsStringAlias.append(1, stateName[0]); 
			}
		}
	PHYCAS_ASSERT(stateListPosAlias.size() == (unsigned)nextStateCode);
	PHYCAS_ASSERT(symbolsStringAlias.size() == (unsigned)nextStateCode);
	this->nativeCMatrix.nObservedStateSets = nextStateCode;
	
	this->nativeCMatrix.nTax = rawMatrix.size();
	this->nativeCMatrix.nChar = (this->nativeCMatrix.nTax == 0 ? 0 : rawMatrix[0].size());
	this->matrixAlias.Initialize(this->nativeCMatrix.nTax, this->nativeCMatrix.nChar);
	nativeCMatrix.matrix = matrixAlias.ptr;
	const unsigned nt = this->nativeCMatrix.nTax;
	const unsigned nc = this->nativeCMatrix.nChar;
	for (unsigned r = 0; r < nt; ++r)
		{
		CIPR_StateSet_t	 * recodedRow = nativeCMatrix.matrix[r];
		const std::vector<int> & rawRowVec = rawMatrix[r];
		//std::cerr <<"row  " << r << '\n';
		PHYCAS_ASSERT(rawRowVec.size() == nc);
		const int * rawRow = &rawRowVec[0];
		for (unsigned c = 0; c < nc; ++c)
			{
			const int rawC = *rawRow++;
			PHYCAS_ASSERT((unsigned)(rawC +  negSCLOffset) < recodeVecLen);
			PHYCAS_ASSERT(rawC >= sclOffset);
			const CIPR_State_t recodedC = recodeArr[rawC];
			PHYCAS_ASSERT(recodedC > -2);
			PHYCAS_ASSERT(recodedC < nextStateCode);
			//std::cerr << "c" << c << ": " << rawC << " => " << (int) recodedC << '\n';
			*recodedRow++ = recodedC;
			}
		}
	nativeCMatrix.symbolsList = symbolsStringAlias.c_str();
	nativeCMatrix.stateListPos = &stateListPosAlias[0];
	nativeCMatrix.stateList = &stateListAlias[0];
	//std::cerr <<"done with DiscreteMatrix ctor\n";

}

/*
void addSymbolForState(string & symbolsForCIPR, const DataStorageType * nclState, const NxsCharactersBlock	& charsBlock, bool verbose)
	{
	try	{
		symbolsForCIPR += datatype.GetStateChar(nclState, false);
		}
	catch (const NxsDataType::XStateNotFound &)
		{
		symbolsForCIPR += ' ';
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



CIPR_Matrix toCIPRMatrix(const ncl::NxsCharactersBlock & inMatrix, const bool verbose)
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

*/
/**
 *	Constructs a native CipresNative::DiscreteMatrix from the native C struct CIPR_Matrix
 *		by deep copy.
 */
DiscreteMatrix::DiscreteMatrix(const CIPR_Matrix & mat)
	:nativeCMatrix(mat),//aliases pointers, but we'll fix this below
	symbolsStringAlias(mat.symbolsList), 
	matrixAlias(mat.nTax, mat.nChar),
	stateListPosAlias(mat.stateListPos, (mat.stateListPos + mat.nObservedStateSets))
	{
	nativeCMatrix.symbolsList = symbolsStringAlias.c_str();
	nativeCMatrix.stateListPos = &stateListPosAlias[0];
	if (mat.nObservedStateSets > 0)
		{
		const unsigned lastStateIndex = nativeCMatrix.stateListPos[nativeCMatrix.nObservedStateSets - 1];
		const unsigned lenAmbigList = lastStateIndex + mat.stateList[lastStateIndex] + 1;
		//	cout << "lenAmbigList = "<< lenAmbigList <<endl;
		stateListAlias.reserve(lenAmbigList);
		cip_copy(mat.stateList, (mat.stateList + lenAmbigList), std::back_inserter(stateListAlias));
		}
	nativeCMatrix.stateList = &stateListAlias[0];
	nativeCMatrix.matrix = matrixAlias.ptr;
	
	// cout << "Matrix in DiscreteMatrix ctor:" << mat.nTax << ' '<< mat.nChar<< endl;
	for (unsigned i = 0; i < mat.nTax; ++i)
		{
		if (mat.nChar > 0)
			cip_copy(mat.matrix[i], mat.matrix[i] + mat.nChar, nativeCMatrix.matrix[i]);
		}
	
	}
			
