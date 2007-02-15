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

//#include "phycas/force_include.h"
#include "phycas/src/oldphycas/const_site_info.hpp"
#include "phycas/src/oldphycas/taxa_manager.hpp"	
#include "phycas/src/ncl/characters/stored_matrix.hpp"	
#include "phycas/src/oldphycas/characters_manager.hpp"
#include "phycas/src/oldphycas/msg_exception.hpp"
using std::vector;
PhoCharactersManager::PhoCharactersManager(PhoTaxaManager & taxaMgr):NxsCharactersManager(taxaMgr){}
/*----------------------------------------------------------------------------------------------------------------------
|	Appends the ordinal codes for the taxon rowIndex. 
*/
void CompressedMatrix::GetTaxonsOrdCodedStates(OrdCodedArrs *oca, unsigned rowIndex) const
	{
	const unsigned offset = rowIndex * dataType.GetNumElementsPerCode();
	VecVecUInt vecOfVecs;
	vecOfVecs.reserve(GetNumPatterns());
	for (PatConstIter pIt = patterns.begin(); pIt != patterns.end(); ++pIt)
		oca->Append(ConvertBitArrToVecUInt(pIt->first + offset, dataType.GetNumElementsPerCode()));
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Fills in a conditional likelihood array with data. Uses 1.0 for observed states, and 0.0 for states that weren't 
|	seen in the row. State array will be copied (to allow for rate heterogeneity). Uses 
|	NxsDataType::FillInCondLikeArray(). The ordering of info is:
|>
|	site0copy0state0 site0copy0state1 ... site0copy1state0 ... site1copy0state0 ...
}<
*/
void CompressedMatrix::FillInTaxonsCondLikeArray(
  double *clArray,						/**< is an array of 1.0 or 0.0 indicating whether or not the state was seen in this row of the matrix */
  unsigned rowInd,						/**< is the index of the row of the matrix (typically the 0-offset taxon number) */
  const unsigned nTimesToCopy) const	/**< is the number of times that each pattern's state array should be copied (e.g. for rate heterogeneity) */
	{
	const unsigned offset = rowInd * dataType.GetNumElementsPerCode();
	const unsigned nStates = dataType.GetNumStates();
	for (PatConstIter pIt = patterns.begin(); pIt != patterns.end(); ++pIt)	
		{
		dataType.FillInCondLikeArray(pIt->first + offset, clArray);
		for (unsigned i = 0; i < nTimesToCopy - 1; ++i)
			{
			nxs_copy_n(clArray, clArray + nStates, nStates);
			clArray += nStates;
			}
		clArray += nStates;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Fills in the ConstSiteInfo structure to describe which patterns can be explained by a invariant model.
|	Uses NxsDataType::FillWithIntersection()
*/
void CompressedMatrix::FillInConstantSiteInfo(
  ConstSiteInfo *csi) const /*on output contains info about which patterns are constant */
	{
	const unsigned nBitCodeUnits = dataType.GetNumElementsPerCode();
	vector<DataPattern> possibleConstantPatterns;
	vector<VecUInt> 		patternIndices; // vector of vectors is horrible storage, but not likely to change size much (there are few types of constant sites)
	DataStorageType    *intersectionCode = new DataStorageType[nBitCodeUnits];	//deleted below (if not "given" to a ConstSiteInfo structure) or in ConstSiteInfo destructor
	PatConstIter patIt = patterns.begin();
	for (unsigned i = 0; patIt != patterns.end(); ++i, ++patIt)
		{
		if (dataType.FillWithIntersection(patIt->first,intersectionCode, nTaxa))
			{
			bool alreadyFound =false;
			for (unsigned j =0; j<possibleConstantPatterns.size(); ++j)
				{
				if (equals_n(intersectionCode, possibleConstantPatterns[j], nBitCodeUnits))
					{
					patternIndices[j].push_back(i);
					alreadyFound = true;
					break;
					}
				}
			if (!alreadyFound)
				{
				VecUInt tempVec = VecUInt(1, i);
				patternIndices.push_back(tempVec);
				possibleConstantPatterns.push_back(intersectionCode);
				intersectionCode = new DataStorageType[nBitCodeUnits];	//deleted below (if not "given" to a ConstSiteInfo structure) or in ConstSiteInfo destructor
				}
			}	
		}
	delete [] intersectionCode;
	csi->nConstStateGrps = (unsigned)possibleConstantPatterns.size();
	NXS_ASSERT(dataType.GetNumStates() < 127); //OrdCodedArrs assumes that a char is enough to hold 
	csi->csStates = OrdCodedArrs(possibleConstantPatterns, nBitCodeUnits, dataType.GetNumStates());
	
	//POL-121803 csi->patternIndices = Compressed2DArray<int>(patternIndices);
	csi->patternIndices = Compressed2DArray<int>(patternIndices, dataType.GetNumStates());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns reference to the DataPatternInfo structure for pattern `i'.
*/
const DataPatternInfo &CompressedMatrix::GetPatternInfo(unsigned i) const
	{
	return patterns[i].second;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Creates and returns a vector of weights for the patterns.  Each pattern's weight will be the sum of the weights of
|	all of the sites that were compressed into the pattern.
*/
vector<double> 	CompressedMatrix::GetPatternWeights() const
	{
	vector<double> w;
	w.reserve(patterns.size());
	for (VecPattern::const_iterator pIt = patterns.begin(); pIt != patterns.end(); ++pIt)
		w.push_back(pIt->second.GetWeight());
	return w;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Deletes data patterns (allocated in PhoCharactersManager::CreateAndInsertDataPatterns)
*/
CompressedMatrix::~CompressedMatrix()
	{
	for (VecPattern::iterator pIt = patterns.begin(); pIt != patterns.end(); ++pIt)
		delete pIt->first;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Uses the PhoCharactersManager::CreateAndInsertDataPatterns to compress and store the characters.
*/
CompressedMatrix::CompressedMatrix(
  const PhoCharactersManager *charMgr, /* character manager that holds the data to be compressed */
  const NxsDataType &desiredType, 	/* datatype that should be used in the encoding/packing */
  const NxsIndexSet &taxaIncluded,	/* indices of the taxa to include */
  const NxsIndexSet &charIndices)	/* indices of the characters to include */
  	:dataType(desiredType),
  	nTaxa(taxaIncluded.size())
  	{
  	DataPatternLess compareF(taxaIncluded.size(), (unsigned char) desiredType.GetNumElementsPerCode()); //@ need to get second arg from the datatype
  	SortedPatterns underConstruction(compareF);
  	charMgr->CreateAndInsertDataPatterns(&underConstruction, desiredType, taxaIncluded, charIndices);
  	patterns.reserve(underConstruction.size());
  	for (SortedPatterns::const_iterator ucIt = underConstruction.begin(); ucIt != underConstruction.end(); ++ucIt)
  		patterns.push_back(PatternAndInfo(ucIt->first, ucIt->second));
  	}
  
	
/*----------------------------------------------------------------------------------------------------------------------
|	Requests a shared pointer to a const ComparessedMatrix
*/
CompressedMatrixShPtr PhoCharactersManager::GetCompressedMatrix(
  const NxsDataType &dt, /*Data type that the characters should be converted to */ //@ need to add functions to assure that conversion is possible
  const NxsIndexSet &taxaIndices) /* indeces of taxa to include in the matrix (all active characters are included) */ 
	{
	NxsIndexSet maskedByActiveTax;
	maskedByActiveTax.SetToIntersection(taxaIndices, NxsTaxaListener::taxaMgr->GetActiveTaxa());
	if (activeMatrixCompressed.first && taxaIndices == activeMatrixCompressed.second)
		return activeMatrixCompressed.first;
	return CompressedMatrixShPtr(new CompressedMatrix(this, dt, taxaIndices, GetActiveCharacters()));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Instantiates arrays to hold pattern descriptions, fills the arrays with patterns for the specified taxa and chars, 
|	and compressses the matrix into a SortedPatterns structure. Called in the CompressedMatrix constructor (which does 
|	most of the work)
*/
void PhoCharactersManager::CreateAndInsertDataPatterns(
  SortedPatterns *dest,					/* on output contains pointers to the patterns codes that must be deleted by the caller */
  const NxsDataType &desiredType,		/* data type to use to store the data */
  const NxsIndexSet &taxaIndices,		/* taxa to be included */
  const NxsIndexSet &charIndices) const /* chars to be included */
  	{
  	if (matrices.empty())
  		return;
  		
  	// Get pointer to first of the matrices stored in the PhoCharactersManager
  	//
  	const ncl::StoredMatrix *currMatrix;
  	VecMatrixShPtr::const_iterator matIt = matrices.begin();
  	currMatrix = (*matIt).get();
  	
  	//@ NEED TO CHECK THAT THE matrices HAVE THE CORRECT TAXA
  	
  	unsigned endCharInThisMat = currMatrix->size() + currMatrix->GetFirstGlobalIndex();
  	
	// cIt goes through the character indices for the current matrix in sequence
  	//
  	NxsIndexSet::const_iterator cIt = charIndices.begin();
  	
  	// Allocate memory for an array to store the pattern for a single character
  	//
  	unsigned lenPatternArray = taxaIndices.size() * desiredType.GetNumElementsPerCode();
  	DataStorageType *nextPattern = new DataStorageType[lenPatternArray];
  			
  	while (cIt != charIndices.end())
  		{
  		while (*cIt >= endCharInThisMat)
  			{
  			++matIt;
  			if (matIt == matrices.end())
  				{
  				delete [] nextPattern; //clean up 
  				NXS_ASSERT(0); //@ we ran out of matrices before getting all of the characters.  Fail silently? exception?
  				return;
  				}
  			currMatrix = (*matIt).get();
  			endCharInThisMat = currMatrix->size() + currMatrix->GetFirstGlobalIndex();
  			}
  			
  		// Copy the pattern represented by this character. The pattern will be represented as bits
  		// e.g. for a DNA site, there are four states (A = 1, C = 2, G = 4, T = 8, ? = 15, - = 0)
  		// 
  		currMatrix->GetPatternFromGlobalIndex(*cIt, nextPattern);
  		
  		SortedPatterns::iterator patIt = dest->find(nextPattern);
  		if (patIt == dest->end())
  			{
  			// This pattern not yet seen
  			(*dest)[nextPattern] = DataPatternInfo(*cIt, 1.0);
  			nextPattern = new DataStorageType[lenPatternArray];
  			}
  		else
  			{
  			// This pattern already known
  			patIt->second.AddIndex(*cIt);
  			}
  		++cIt;
  		}
  		
  	delete [] nextPattern;
  	}
