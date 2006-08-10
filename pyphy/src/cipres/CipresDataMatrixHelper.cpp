#include "pyphy/src/cipres/CipresDataMatrixHelper.h"
#include "pyphy/src/cipres/util_copy.hpp"
using CipresNative::DiscreteMatrix;
using std::string;
using std::vector;
using std::cout;
using std::endl;

#define MISSING_INCLUDES_GAP	/* If not defined, then missing does not include the "gap" state (-1).
								   Ordinarily, this SHOULD be defined */

#if !defined(NO_IDL_TYPES)
CIPR_Datatypes CipresNative::corbaToNativeDatatype(const CipresIDL_api1::DiscreteDatatypes dt)
	{
	switch (dt)
		{
		case CipresIDL_api1::DNA_DATATYPE:		return CIPR_DNA_Datatype;
		case CipresIDL_api1::RNA_DATATYPE:		return CIPR_RNA_Datatype;
		case CipresIDL_api1::AA_DATATYPE:		return CIPR_AA_Datatype;
		case CipresIDL_api1::CODON_DATATYPE:		return CIPR_Codon_Datatype;
		case CipresIDL_api1::GENERIC_DATATYPE:	return CIPR_Generic_Datatype;
		default:
			assert(0);
		}
		//unreachable
	return CIPR_Generic_Datatype;
	}

CipresIDL_api1::DiscreteDatatypes CipresNative::nativeToCorbaDatatype(const CIPR_Datatypes dt)
	{
	switch (dt)
		{
		case CIPR_DNA_Datatype: 	return CipresIDL_api1::DNA_DATATYPE;
		case CIPR_RNA_Datatype: 	return CipresIDL_api1::RNA_DATATYPE;
		case CIPR_AA_Datatype: 		return CipresIDL_api1::AA_DATATYPE;
		case CIPR_Codon_Datatype:	return CipresIDL_api1::CODON_DATATYPE;
		case CIPR_Generic_Datatype:	return CipresIDL_api1::GENERIC_DATATYPE;
		default:
			assert(0);
		}
		//unreachable
	return CipresIDL_api1::GENERIC_DATATYPE;
	}

/**
 *	Creates a CipresIDL_api1::DataMatrix form of the matrix suitable for return as a CORBA return type 
 *		(just like the _retn function of a *_var class).
 */
CipresIDL_api1::DataMatrix *DiscreteMatrix::_retn() const
	{
	CipresIDL_api1::DataMatrix_var corbaMatrix(new CipresIDL_api1::DataMatrix);
	corbaMatrix->m_symbols = CORBA::string_dup(nativeCMatrix.symbolsList);
	corbaMatrix->m_numStates = nativeCMatrix.nStates;
	corbaMatrix->m_numCharacters = nativeCMatrix.nChar;
	corbaMatrix->m_datatype = nativeToCorbaDatatype(nativeCMatrix.datatype);
	corbaMatrix->m_matrix.length(nativeCMatrix.nTax);
	for (unsigned i = 0; i < nativeCMatrix.nTax; ++i)
		{
		CipresIDL_api1::Characters & row = corbaMatrix->m_matrix[i];
		const CIPR_StateSet_t * sourcePtr = nativeCMatrix.matrix[i];
		row.length(nativeCMatrix.nChar);
		for (unsigned j = 0; j < nativeCMatrix.nChar; ++j)
			row[j] = CORBA::Short(*sourcePtr++);
		}
	corbaMatrix->m_charStateLookup.length(nativeCMatrix.nObservedStateSets);
	//	cout << " nativeCMatrix.nObservedStateSets = "<< nativeCMatrix.nObservedStateSets<<endl;
	for (unsigned k = 0; k < nativeCMatrix.nObservedStateSets; ++k)
		{
		CipresIDL_api1::StateSet & stateSet = corbaMatrix->m_charStateLookup[k];
		const unsigned stateListInd = nativeCMatrix.stateListPos[k];
		const unsigned nStatesInSet = nativeCMatrix.stateList[stateListInd];
		stateSet.length(nStatesInSet);
		//	cout << "State set " << k << endl;
		for (unsigned m = 0; m < nStatesInSet; ++m)
			{
			stateSet[m] = nativeCMatrix.stateList[stateListInd + 1 + m];
			//	cout << " stateSet["<< m <<"] = "<< stateSet[m];
			}
		//	cout <<endl;
		}
	return corbaMatrix._retn();
	}

/**
 *	Constructs a native CipresNative::DiscreteMatrix from a corba type (CipresIDL_api1::DataMatrix)
 *	Throws CORBA exceptions on errors (CORBA::NO_MEMORY and CipresIDL_api1::CharStateLookup)
 */
DiscreteMatrix::DiscreteMatrix(const CipresIDL_api1::DataMatrix & corbaMatrix)
	{
	try	
		{
		nativeCMatrix.datatype = corbaToNativeDatatype(corbaMatrix.m_datatype);
		nativeCMatrix.nStates = corbaMatrix.m_numStates;
		nativeCMatrix.nChar = corbaMatrix.m_numCharacters;
		nativeCMatrix.nTax = corbaMatrix.m_matrix.length();
		matrixAlias.Initialize(nativeCMatrix.nTax, nativeCMatrix.nChar);
		nativeCMatrix.matrix = matrixAlias.ptr;
		for (unsigned i = 0; i < nativeCMatrix.nTax; ++i)
			{
			const CipresIDL_api1::Characters & x_i = corbaMatrix.m_matrix[i];
			for (unsigned j = 0; j < nativeCMatrix.nChar; j++)
				nativeCMatrix.matrix[i][j] = static_cast<CIPR_StateSet_t>(x_i[j]);
			}
		symbolsStringAlias = string(corbaMatrix.m_symbols);
		nativeCMatrix.symbolsList = symbolsStringAlias.c_str();
		
		const CipresIDL_api1::CharStateLookup stateLookup = corbaMatrix.m_charStateLookup;
		const unsigned nObservedStateSets = stateLookup.length();
		stateListPosAlias.reserve(nObservedStateSets);
		nativeCMatrix.nObservedStateSets = nObservedStateSets;
		stateListPosAlias.reserve(nObservedStateSets);
		if (nObservedStateSets <= nativeCMatrix.nStates)
			throw CipresIDL_api1::BadArgException("DataMatrix Error: m_charStateLookup.length() <= m_numStates");
		unsigned stateSetIndex = 0;
		for (; stateSetIndex < nativeCMatrix.nStates; ++stateSetIndex)
			{
			stateListAlias.push_back(1);
			stateListAlias.push_back(static_cast<CIPR_State_t>(stateSetIndex));
			stateListPosAlias.push_back(2*stateSetIndex);
			if (stateSetIndex != (unsigned) stateLookup[stateSetIndex][0])
				throw CipresIDL_api1::BadArgException("DataMatrix Error: m_charStateLookup[i][0] != i for all i <= m_numStates");
			}
		stateListPosAlias.push_back(2*stateSetIndex);
#		if defined(MISSING_INCLUDES_GAP)
			stateListAlias.push_back(nativeCMatrix.nStates + 1);
			stateListAlias.push_back(-1);
#		else
			stateListAlias.push_back(nativeCMatrix.nStates);
#		endif
		assert(nativeCMatrix.nStates < 127); //need to fit into int8_t
		for (CIPR_State_t i = 0; i < static_cast<CIPR_State_t>(nativeCMatrix.nStates); ++i)
			stateListAlias.push_back(i);
		++stateSetIndex;
		for (; stateSetIndex < nObservedStateSets;  ++stateSetIndex)
			{
			stateListPosAlias.push_back(stateListAlias.size());
			
			const CipresIDL_api1::StateSet & currStateSet = stateLookup[stateSetIndex];
			CIPR_State_t currNStatesObserved = static_cast<CIPR_State_t>(currStateSet.length());
			if (currNStatesObserved == 1)
				throw CipresIDL_api1::BadArgException("DataMatrix Error: unambiguous state set in m_charStateLookup");
			if (currNStatesObserved > static_cast<CIPR_State_t>(nativeCMatrix.nStates))
				throw CipresIDL_api1::BadArgException("DataMatrix Error: list of ambiguous states is longer than m_numStates");
			
			stateListAlias.push_back(currNStatesObserved);
			for (CIPR_State_t i = 0; i < currNStatesObserved; ++i)
				{
				const int currState = currStateSet[i];
				if (currState < 0)
					throw CipresIDL_api1::BadArgException("DataMatrix Error: negative number found in m_charStateLookup (currently, CIPR_Matrix does not support ambiguity with gaps codes).");
				if (currState >= (int)nativeCMatrix.nStates)
					throw CipresIDL_api1::BadArgException("DataMatrix Error: number >= m_numStates found in m_charStateLookup");
				stateListAlias.push_back(static_cast<CIPR_State_t>(currState));
				}
			}
		nativeCMatrix.stateList = &stateListAlias[0];
		nativeCMatrix.stateListPos = &stateListPosAlias[0];
		}
	catch (std::bad_alloc &)
		{
		throw CORBA::NO_MEMORY();
		}
	}
#endif

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
	/*
	cout << "DiscreteMatrix ctor stateListPosAlias" << endl;
	for (unsigned i = 0; i < stateListPosAlias.size(); ++i)
		cout << '\t'<< int(stateListPosAlias[i]) << endl;
	cout << "mat.nObservedStateSets = "<< mat.nObservedStateSets << ' '<< nativeCMatrix.nObservedStateSets <<endl;
	*/
	if (mat.nObservedStateSets > 0)
		{
		const unsigned lastStateIndex = nativeCMatrix.stateListPos[nativeCMatrix.nObservedStateSets - 1];
		const unsigned lenAmbigList = lastStateIndex + mat.stateList[lastStateIndex] + 1;
		//	cout << "lenAmbigList = "<< lenAmbigList <<endl;
		stateListAlias.reserve(lenAmbigList);
		cip_copy(mat.stateList, (mat.stateList + lenAmbigList), std::back_inserter(stateListAlias));
		}
	/*
	cout << "DiscreteMatrix ctor stateListAlias" << endl;
	for (unsigned i = 0; i < stateListAlias.size(); ++i)
		cout << '\t'<< int(stateListAlias[i]) << endl;
	*/
	nativeCMatrix.stateList = &stateListAlias[0];
	nativeCMatrix.matrix = matrixAlias.ptr;
	
	// cout << "Matrix in DiscreteMatrix ctor:" << mat.nTax << ' '<< mat.nChar<< endl;
	for (unsigned i = 0; i < mat.nTax; ++i)
		{
		if (mat.nChar > 0)
			cip_copy(mat.matrix[i], mat.matrix[i] + mat.nChar, nativeCMatrix.matrix[i]);
		/*
		cout << i << '\t';
		for (unsigned j = 0; j < mat.nChar; ++j)
			cout <<  int(nativeCMatrix.matrix[i][j]);
		cout << endl;
		*/
		}
	
	}
			
