#ifndef NCL_NXSDATATYPE_INLINES_H
#define NCL_NXSDATATYPE_INLINES_H

#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/misc/nxs_data_type.hpp"
#include "pyphy/src/ncl/misc/nxs_discrete_matrix.hpp"
#include "pyphy/src/ncl/nxs_token.hpp"



	
/*--------------------------------------------------------------------------------------------------------------------------
|	detemines whether the character being read is a {, (, or a character symbol, and calls the appropriate function 
|	to continue parsing.
*/
template<class T> inline void NxsDataType::ReadNextState(
  T *dataMatrix,
  NxsToken& token, 
  unsigned i, 	/*taxon index */
  unsigned j) const  /*character index, pass UINT_MAX to skip this character (for eliminated characters) */
	{
	// this should be the state for taxon i and character j
	//
	const char &ch(token.GetTokenReference()[0]);
	if (ch == '(')
		ReadNextMultipleStates<T>(dataMatrix, i, j, token, true);
	else if (ch =='{')
		ReadNextMultipleStates<T>(dataMatrix, i, j, token, false);
	else if (j!= UINT_MAX)	
		dataMatrix->SetStateChar(i, j, ch, token);
	}

// don't call with eliminated chars
// too long to inline, but called by very few functions in a speed critical portion of code, so we'll go ahead and bloat a bit
inline void NxsDataType::AddStateChar(
  ncl::NxsDiscreteMatrix * dataMatrix,
  unsigned i, 
  unsigned j, 
  char ch,
  const NxsToken &token) const
  	{
  	NXS_ASSERT(j != UINT_MAX && i != UINT_MAX);
  	const DataStorageType *codedState = GetStateInBitCode(ch);
	if (codedState != NULL)
  		dataMatrix->AddState(i, j, codedState);
  	else if (ch == missingSymbol) 
		dataMatrix->SetMissing(i, j);	// set missing and adding the missing code should be equivalent
	else if(ch == gapSymbol) 
		dataMatrix->AddGap(i, j); 
	else if (ch == matchSymbol)
		dataMatrix->AddState(i, j, dataMatrix->GetState(0, j));
	else
		{
		std::map<char, EquateValInfo>::const_iterator eIt = equateCodes.find(ch);
  		if (eIt == equateCodes.end())
  			{
  			std::string errormsg;
  			errormsg << "State specified (" << token.GetTokenReference() << ") for taxon " << (i+1) << ", character " << (j+1) << ", not found in list of valid symbols";
			throw NxsException(errormsg, token);
			}
		else
			{
			dataMatrix->AddState(i, j, eIt->second.first);
			if (eIt->second.second == 1)
				dataMatrix->SetPolymorphic(i,j);
			}
		}
	} 
//don't call with eliminated characters
// too long to inline, but called by very few functions in a speed critical portion of code, so we'll go ahead and bloat a bit
inline void NxsDataType::AddStateChar(
  ncl::NxsUnalignedMatrix *dataMatrix,
  unsigned i, 
  unsigned j, 
  char ch,
  const NxsToken &token) const
  	{
  	PHYCAS_ASSERT(j != UINT_MAX && i != UINT_MAX);
  	const DataStorageType *codedState = GetStateInBitCode(ch);
  	if (codedState != NULL)
  		dataMatrix->AddState(i, j, codedState);
  	else if (ch == missingSymbol) 
		dataMatrix->SetMissing(i, j);	// set missing and adding the missing code should be equivalent
	else
		{
		std::map<char, EquateValInfo>::const_iterator eIt = equateCodes.find(ch);
  		if (eIt == equateCodes.end())
  			{
  			std::string errormsg;
  			errormsg << "State specified (" << token.GetTokenReference() << ") for taxon " << (i+1) << ", character " << (j+1) << ", not found in list of valid symbols";
			throw NxsException(errormsg, token);
			}
		else
			{
			dataMatrix->AddState(i, j, eIt->second.first);
			if (eIt->second.second == 1)
				dataMatrix->SetPolymorphic(i,j);
			}
		}
	} 
// too long to inline, but called by very few functions in a speed critical portion of code, so we'll go ahead and bloat a bit
inline void ncl::NxsDiscreteMatrix::SetStateChar(
  unsigned i, 
  unsigned j, 
  char ch,
  const NxsToken &token)
  	{
  	PHYCAS_ASSERT(j != UINT_MAX && i != UINT_MAX);
  	const DataStorageType *codedState = dataType.GetStateInBitCode(ch);
  	if (codedState != NULL)
  		SetState(i, j, codedState);
  	else if (ch == dataType.missingSymbol)
		SetMissing(i, j);
	else if(ch == dataType.gapSymbol) 
		SetGap(i, j);
	else if (ch == dataType.matchSymbol)
		CopyStatesFromFirstRow(i, j);
	else
		{
		std::map<char, NxsDataType::EquateValInfo>::const_iterator eIt = dataType.equateCodes.find(ch);
		if (eIt == dataType.equateCodes.end())
			{
			std::string errormsg;
			errormsg << "State specified (" << token.GetTokenReference() << ") for taxon " << (i+1) << ", character " << (j+1) << ", not found in list of valid symbols";
			throw NxsException(errormsg, token);
			}
		else
			{
			SetState(i, j, eIt->second.first);
			if (eIt->second.second == 1)
				SetPolymorphic(i,j);
			}
		}
	} 

// too long to inline, but called by very few functions in a speed critical portion of code, so we'll go ahead and bloat a bit
inline void ncl::NxsUnalignedMatrix::SetStateChar(
  unsigned i, 
  unsigned j, 
  char ch,
  const NxsToken &token)
  	{
  	PHYCAS_ASSERT(j != UINT_MAX && i != UINT_MAX);
  	const DataStorageType *codedState =dataType. GetStateInBitCode(ch);
  	if (codedState != NULL)
  		SetState(i, j, codedState);
  	else if (ch == dataType.missingSymbol) 
		SetMissing(i, j);
	else
		{
		std::map<char, NxsDataType::EquateValInfo>::const_iterator eIt = dataType.equateCodes.find(ch);
  		if (eIt == dataType.equateCodes.end())
  			{
  			std::string errormsg;
  			errormsg << "State specified (" << token.GetTokenReference() << ") for taxon " << (i+1) << ", character " << (j+1) << ", not found in list of valid symbols";
			throw NxsException(errormsg, token);
			}
		else
			{
			SetState(i, j, eIt->second.first);
			if (eIt->second.second == 1)
				SetPolymorphic(i,j);
			}
		}
	} 
template <class T> void BitShiftAndUnion(
  T* arrStart,
  unsigned lenOfArray)
  {
  const T kMaxBitOfTypeT = (T) (((T) ~ 0) - (((T) ~ 0) >> 1));
  bool carryTheOne = false;
  for (unsigned i = 0; i < lenOfArray; ++i)
  	{
  	if (carryTheOne)
  		arrStart[i] |= (T) 1;
  	carryTheOne = ((arrStart[i] & kMaxBitOfTypeT) != 0);
  	arrStart[i] |= (arrStart[i] << 1);
  	}
  }

template <class T> void BitShiftAndUnion_n(
  T* arrStart,
  unsigned lenOfArray,
  unsigned nTimesToRepeat)
  {
  for (unsigned i = 0; i < nTimesToRepeat; ++i)
  	BitShiftAndUnion<T>(arrStart, lenOfArray);
  }

template <class T> void NxsDataType::AddStateRange(
  T *dataMatrix,
  unsigned i, 
  unsigned j, 
  char fromC,
  char toC,
  const NxsToken &token) const
  	{
  	const DataStorageType *codedFrom = GetStateInBitCode(fromC);
  	const DataStorageType *codedTo = GetStateInBitCode(toC);
  	
  	if (*codedFrom != MAX_NSTATES_ALLOWED && *codedTo != MAX_NSTATES_ALLOWED)
		{
		unsigned ordFrom = IndexOfFirstOnBit(codedFrom, nCodeUnitsPerCharacter) ;
  		unsigned ordTo = IndexOfFirstOnBit(codedTo, nCodeUnitsPerCharacter) ;
  		PHYCAS_ASSERT(ordFrom != UINT_MAX && ordTo != UINT_MAX);
  		if (ordFrom < ordTo)
			{
  			if (j != UINT_MAX)
				{
				// fill in scratchSpace with all codes between codedFrom and codedTo;
				//
				nxs_copy(codedFrom, codedFrom + nCodeUnitsPerCharacter, scratchSpace);
				BitShiftAndUnion_n<DataStorageType>(scratchSpace, nCodeUnitsPerCharacter, ordTo - ordFrom);
				dataMatrix->AddState(i, j, scratchSpace);
				}
			return;
			}
		}
	std::string errormsg;
	errormsg << "Invalid state range " << fromC << " ~ " << toC;
	throw NxsException(errormsg, token);
	}

inline std::string ncl::NxsDiscreteMatrixBase::GetStatesAsNexusString(unsigned i, unsigned j, bool expandEquates) const
	{
	return dataType.GetStatesAsNexusString(GetState(i, j), IsPolymorphic(i, j), expandEquates);
	}


#endif

