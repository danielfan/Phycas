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
#include "phypy/src/ncl/misc/nxs_discrete_matrix.hpp"
#include "phypy/src/ncl/misc/nxs_data_type_inl.hpp"
#include <cmath>

using ncl::NxsDiscreteMatrixBase;
using ncl::NxsUnalignedMatrix;
using ncl::NxsDiscreteMatrix;

NxsDiscreteMatrixBase::NxsDiscreteMatrixBase(
  NxsDataType DatType)
  	:dataType(DatType)
  	{
#	if defined(C_FUNCS_IN_STD_NAMESPACE)
		using std::floor;
#	endif
  	unsigned nBitsInStoreType = 8*sizeof(DataStorageType);
  	nUnitsPerElement = 1 + (unsigned) floor (((double) (dataType.GetNumStates() - 1) )/ (double) nBitsInStoreType);
  	}
 
NxsUnalignedMatrix::NxsUnalignedMatrix(
  unsigned rows, 
  NxsDataType DatType)
  	:NxsDiscreteMatrixBase(DatType)
  	{
  	DiscreteDataRow dr;
  	dataMat.reserve(rows);
  	for (unsigned i = 0; i < rows; ++i)
  		dataMat.push_back(dr);
  	}
  
void NxsDiscreteMatrixBase::RemoveRow(
  unsigned rowN)
  	{
  	if (rowN < dataMat.size())
  		{
  		DiscreteDataMatrix::iterator matIt = dataMat.begin();
  		for (; rowN > 0; --rowN)
  			++matIt;
  		dataMat.erase(matIt);
  		}
  	}

//	using the swap trick to get rid of extra rows
//
void NxsDiscreteMatrixBase::TrimExcessCapacity()
  	{
  	DiscreteDataMatrix(dataMat).swap(dataMat);
  	}


NxsDiscreteMatrix::NxsDiscreteMatrix(
  unsigned rows, 
  unsigned cols, 
  NxsDataType DatType)
  	:NxsDiscreteMatrixBase(DatType),
  	nCols(cols)
  	{
  	DiscreteDataRow dr;
  	dr.reserve(nCols * nUnitsPerElement);;
  	const DataStorageType  *tempMissingCode = dataType.GetMissingBitCode();
  	for (unsigned j = 0;j < nCols; ++j)
  		{
  		for (unsigned k = 0; k < nUnitsPerElement; ++k)
  			dr.push_back(tempMissingCode[k]);
  		}
  	dataMat.reserve(rows);
  	for (unsigned i = 0; i < rows; ++i)
  		dataMat.push_back(dr);
  	}
 

/*--------------------------------------------------------------------------------------------------
|	will always be zero for single DataStorageTypes
*/
inline unsigned	NxsDiscreteMatrixBase::UnitToAddTo(const unsigned &ordination) const
	{
	const unsigned nBitsInStoreType = 8*sizeof(DataStorageType);
  	NXS_ASSERT (ordination / nBitsInStoreType <= nUnitsPerElement);
  	return ordination / nBitsInStoreType;
	}
	
inline unsigned UnsignedBinaryPatternFromOrd(unsigned ordination)
	{
	const unsigned nBitsInStoreType = 8*sizeof(DataStorageType);
	while (ordination >= nBitsInStoreType)
		ordination -= nBitsInStoreType;
	unsigned retValue = 1U;
	return (retValue << ordination);
  	}
	
void NxsDiscreteMatrixBase::AddStateIndex(
  unsigned rowN, 
  unsigned colN,
  unsigned ordinationOfState)
	{
	NXS_ASSERT( IsValidPutElement(rowN, colN));
  	unsigned unitToAddTo = UnitToAddTo(ordinationOfState);
	unsigned binaryPatternToAdd = UnsignedBinaryPatternFromOrd(ordinationOfState);
  	
  	dataMat[rowN][nUnitsPerElement*colN + unitToAddTo] |= binaryPatternToAdd;
	}
	
/**
 * @method Flush [void:public]
 *
 * Deletes all cells of dataMat and resets nrows and ncols to 0.
 */
void NxsDiscreteMatrixBase::Flush()
	{
	dataMat.clear();
	}
	
void NxsDiscreteMatrix::Flush()
	{
	NxsDiscreteMatrixBase::Flush();
	nCols = 0;
	}
	

void NxsDiscreteMatrix::DebugSaveMatrix(
  std::ostream& out) const
	{
	out << "\nnrows = " << GetNumRows() << "\nncols = " << GetNumColumns();
	for(unsigned i = 0; i < GetNumRows(); ++i) 
		{
		out << "\n";
		for(unsigned j = 0; j < GetNumColumns(); ++j) 
			{
			if(IsMissing(i, j))
				{
				for ( unsigned k = 0; k < nUnitsPerElement; ++k)
				   out << '?';
				}
			else if(IsGap(i, j))
				{
				for ( unsigned k = 0; k < nUnitsPerElement; ++k)
				   out << '-';
				}
			else
				{
				for ( unsigned k = 0; k < nUnitsPerElement; ++k)
				   out << GetState(i, j)[k];
				}
			}
		}
	out << std::endl;
	}

bool NxsUnalignedMatrix::PadRowZeroLast(
  unsigned rowN, 
  unsigned colN)
  	{
  	NXS_ASSERT(rowN < GetNumRows());
  	const DataStorageType  *tempMissingCode = dataType.GetMissingBitCode();
  	unsigned nToAdd = colN + 1 - (dataMat[rowN].size()/ nUnitsPerElement);
  	for  (; nToAdd > 1; --nToAdd)
		{
		for (unsigned i = 0; i < nUnitsPerElement; ++i)
			dataMat[rowN].push_back(tempMissingCode[i]);
		}
	if (nToAdd > 0)
		{
		for (unsigned i = 0; i < nUnitsPerElement; ++i)
			dataMat[rowN].push_back(0);
		return true;
		}
	return false;
	}


