#ifndef NCL_NXSDISCRETEMATRIX_H
#define NCL_NXSDISCRETEMATRIX_H

#include "pyphy/src/ncl/misc/nxs_data_type.hpp"
#include "pyphy/src/ncl/misc/algorithm_extensions.hpp"

typedef std::vector<DataStorageType> DiscreteDataRow;
typedef std::vector<DiscreteDataRow> DiscreteDataMatrix;

typedef std::set<unsigned> IndSet;
typedef std::map<unsigned, IndSet> TwoDLookUp;

namespace ncl
{
/*--------------------------------------------------------------------------------------------------
|	Stores a matrix of data.  NxsDiscreteMatrixBase can handle an arbitrary number of states.  
|	Data is encoded with bits standing for states.  A dataType object is needed to convert the codes
|	to characters (get info on which bit stands for which state).
|	The fundamental storage type, DataStorageType, may be too small to fit every type of character.
|	If this is the case multiple adjacent DataStorageTypes will be used, so each state is stored in
|	an array.
|	If access speed is critical GetState(taxonIndex, 0) can be called to get a const pointer to the 
|	beginning of the taxon's character array.
*/
class NxsDiscreteMatrixBase
	{
	public:
		NxsDiscreteMatrixBase(NxsDataType DatType);
		virtual ~NxsDiscreteMatrixBase() {}
		
		
		//	Accessors
		//
		const NxsDataType 	   &GetDataType() const;
		unsigned 				GetNumRows()	const;
		const DataStorageType  *GetState(unsigned rowN, unsigned colN) const;
		void					GetPattern(unsigned colN, DataStorageType *p) const;
		virtual bool			IsMissing(unsigned rowN, unsigned colN) const;
		bool					IsPolymorphic(unsigned rowN, unsigned colN) const;
		std::string				GetStatesAsNexusString(unsigned i, unsigned j, bool expandEquates) const;
		 
		//  Modifiers
		//
		virtual void 			AddStateIndex(unsigned rowN, unsigned colN,unsigned);
		virtual void 			AddState(unsigned rowN, unsigned colN,const DataStorageType *value);
		virtual void	 		Flush();
		void 					SetStateIndex(unsigned rowN, unsigned colN,unsigned);
		virtual void 			SetState(unsigned rowN, unsigned colN,const DataStorageType *value);
		virtual void			SetMissing(unsigned rowN, unsigned colN);
		void					SetPolymorphic(unsigned rowN, unsigned colN);
		void					RemoveRow(unsigned);
		void					TrimExcessCapacity();
		virtual void			ZeroState(unsigned rowN, unsigned colN);
		virtual void	 		SetStateChar(unsigned rowN, unsigned colN, char ch, const NxsToken &) = 0;
		protected:
	
		virtual bool	 		IsValidGetElement(unsigned rowN, unsigned colN) const = 0;
		virtual bool	 		IsValidPutElement(unsigned rowN, unsigned colN) const = 0;
		unsigned 				UnitToAddTo(const unsigned &ordination) const;
		
		typedef std::set<unsigned> IndSet;
		typedef std::map<unsigned, IndSet> TwoDLookUp;
		
		NxsDataType   		dataType;
		unsigned			nUnitsPerElement;	/* this is stride through the character bit code array.  It will be the same as NxsDataType's nCodeUnitsPerCharacter, but is duplicated because it is used frequently here*/
		DiscreteDataMatrix  dataMat;
		TwoDLookUp			polymorphic;	
	};

/*--------------------------------------------------------------------------------------------------
|	Stores aligned datamatrices.
|	Note that gapped characters have all zero bits, missing data has all bits of the datatype set to
|	1
*/
class NxsDiscreteMatrix : public NxsDiscreteMatrixBase
	{
	public:
		NxsDiscreteMatrix(unsigned rows, unsigned cols, NxsDataType DatType);
		
		void 					DebugSaveMatrix(std::ostream& out) const;
		
		//	Accessors
		//
		unsigned 				GetNumColumns() const;
		bool					IsGap(unsigned rowN, unsigned colN) const;
		const IndSet		  & GetGapMatrixRow(unsigned rowN) const;
		virtual bool			IsMissing(unsigned rowN, unsigned colN) const;
		 
		//  Modifiers
		//
		void 					AddGap(unsigned rowN, unsigned colN);
		void 					CopyStatesFromFirstRow(unsigned rowN, unsigned colN);
		void	 				Flush();
		void					SetGap(unsigned rowN, unsigned colN);
		virtual void			SetMissing(unsigned rowN, unsigned colN);
		
		void	 				SetStateChar(unsigned rowN, unsigned colN, char ch, const NxsToken &);
	protected:
	
		
		bool	 				IsValidGetElement(unsigned rowN, unsigned colN) const;
		bool	 				IsValidPutElement(unsigned rowN, unsigned colN) const;
		
		unsigned 			nCols;
		TwoDLookUp			gapped;
		const IndSet		emptyIndSet; // for return from GetGapMatrixRow when there are no gaps in a row
	};
	
		
/*--------------------------------------------------------------------------------------------------
|	Stores unaligned datamatrices.  Note GetNumColumns takes an argument (the row number) because
|	# of columns varies.  
*/
class NxsUnalignedMatrix : public NxsDiscreteMatrixBase
	{
	public:
		NxsUnalignedMatrix(unsigned rows, NxsDataType DatType);
		
		void 					DebugSaveMatrix(std::ostream& out) const;
		
		//	Accessors
		//
		unsigned 				GetNumColumns(unsigned rowN) const;
		 
		//  Modifiers
		//
		void 					AddStateIndex(unsigned rowN, unsigned colN,unsigned);
		void 					AddState(unsigned rowN, unsigned colN,const DataStorageType *value);
		void 					SetState(unsigned rowN, unsigned colN,const DataStorageType *value);
		void					ZeroState(unsigned rowN, unsigned colN);

		void	 				SetStateChar(unsigned rowN, unsigned colN, char ch, const NxsToken &);
	
	protected:
		bool	 				IsValidGetElement(unsigned rowN, unsigned colN) const;
		bool	 				IsValidPutElement(unsigned rowN, unsigned colN) const;
		bool					PadRowZeroLast(unsigned rowN, unsigned colN);
		
	};
	
	
typedef boost::shared_ptr<NxsDiscreteMatrixBase> DiscreteMatrixBaseShPtr;
typedef boost::shared_ptr<NxsDiscreteMatrix> DiscreteMatrixShPtr;
typedef boost::shared_ptr<NxsUnalignedMatrix> UnalignedMatrixShPtr;
		
/*--------------------------------------------------------------------------------------------------
|	returns a const reference to the NxsDataType object that can be used to interpret the bit codes
|	that the Matrix stores
*/
inline const NxsDataType &NxsDiscreteMatrixBase::GetDataType() const
	{
	return dataType;
	}
		
/*--------------------------------------------------------------------------------------------------
|	returns a number of rows (typically the number of taxa)
*/
inline unsigned NxsDiscreteMatrixBase::GetNumRows()	const 
	{
	return (unsigned)dataMat.size();
	}
	
/*--------------------------------------------------------------------------------------------------
|	returns a number of columns (typically the number of characters)
*/
inline unsigned NxsDiscreteMatrix::GetNumColumns() const
	{
	return nCols;
	}

/*--------------------------------------------------------------------------------------------------
|	returns a number of columns for `rowN'
*/
inline unsigned NxsUnalignedMatrix::GetNumColumns(unsigned rowN) const
	{
	assert(rowN < dataMat.size());
	return (unsigned)dataMat[rowN].size()/nUnitsPerElement;
	}

/*--------------------------------------------------------------------------------------------------
|	returns true if the args describe a valid cell in the matrix to be read from.  //@temp
|	we're currently not keeping track of which elements have been initialized.  So this function is 
|	bogus (but we are reading in the entire matrix in characters blocks so this function behaves
|	correctly in that context)
*/
inline bool NxsDiscreteMatrix::IsValidGetElement(
  unsigned rowN,
  unsigned colN) const
  	{
  	return (rowN < GetNumRows() && colN < GetNumColumns());
  	}
  	
/*--------------------------------------------------------------------------------------------------
|	returns true if the args describe a valid cell in the matrix whose state can be set.
*/
inline bool NxsDiscreteMatrix::IsValidPutElement(
  unsigned rowN,
  unsigned colN) const
  	{
  	return IsValidGetElement(rowN, colN);
  	}
  	
/*--------------------------------------------------------------------------------------------------
|	returns true if the args describe a valid cell in the matrix to be read from.  //@temp
|	we're currently not keeping track of which elements have been initialized.  So this function is 
|	bogus (but we are reading in the entire matrix in unaligned blocks so this function behaves
|	correctly in that context)
*/
inline bool NxsUnalignedMatrix::IsValidGetElement(
  unsigned rowN,
  unsigned colN) const
  	{
  	return (rowN < GetNumRows() && colN < GetNumColumns(rowN));
  	}
	
	
/*--------------------------------------------------------------------------------------------------
|	copies the bit code from inValue into cell [rowN][colN]
*/
inline void NxsDiscreteMatrixBase::SetState(
  unsigned rowN, 
  unsigned colN, 
  const DataStorageType *inValue)
  	{
  	assert(IsValidPutElement(rowN, colN));
  	nxs_copy(inValue, (inValue + nUnitsPerElement), &dataMat[rowN][nUnitsPerElement*colN]);
  	}
  	
/*--------------------------------------------------------------------------------------------------
|	copies the bit code from inValue into cell [rowN][colN]
*/
inline void NxsUnalignedMatrix::SetState(
  unsigned rowN, 
  unsigned colN,
  const DataStorageType *inValue)	
	{
	PadRowZeroLast(rowN, colN);
	NxsDiscreteMatrixBase::SetState(rowN, colN, inValue);
  	}
  	
/*--------------------------------------------------------------------------------------------------
|	adds the bits in the bit code from inValue into cell [rowN][colN]
*/
inline void NxsDiscreteMatrixBase::AddState(
  unsigned rowN, 
  unsigned colN, 
  const DataStorageType *inValue)
  	{
  	assert(IsValidPutElement(rowN, colN));
  	for (unsigned i = 0; i < nUnitsPerElement; ++i)
  		dataMat[rowN][nUnitsPerElement*colN + i] |= inValue[i];
  	}
  	
/*--------------------------------------------------------------------------------------------------
|	sets the bits of the cell [rowN][colN] to 0
*/
inline void NxsDiscreteMatrixBase::ZeroState(
  unsigned rowN, 
  unsigned colN)
  	{
  	assert(IsValidPutElement(rowN, colN));
  	for (unsigned i = 0; i < nUnitsPerElement; ++i)
  		dataMat[rowN][nUnitsPerElement*colN + i] = (DataStorageType) 0;
  	}

/*--------------------------------------------------------------------------------------------------
|	same as SetState, but takes the index of the state as opposed to the bitcode
*/
inline void NxsDiscreteMatrixBase::SetStateIndex(unsigned rowN, unsigned colN,unsigned indexOfState)
	{
	ZeroState(rowN, colN);
	AddStateIndex(rowN, colN, indexOfState);
	}
  	
/*--------------------------------------------------------------------------------------------------
|	Sets the state for cell[rowN][colN] to the state in cell[0][colN]
*/
inline void NxsDiscreteMatrix::CopyStatesFromFirstRow(
  unsigned rowN, 
  unsigned colN)
  	{
  	assert(IsValidPutElement(rowN, colN));
  	if (IsGap(0,colN))
  		SetGap(rowN,colN);
  	else
  		SetState(rowN, colN, GetState(0, colN));
  	}

/*--------------------------------------------------------------------------------------------------
|	Returns the internal pointer to the bitcode for cell[rowN][colN]
*/
inline const DataStorageType *NxsDiscreteMatrixBase::GetState(
  unsigned rowN, 
  unsigned colN) const
  	{
  	assert(IsValidGetElement(rowN, colN));
  	return ((&dataMat[rowN][0]) + (colN * nUnitsPerElement));
  	}
  	
/*--------------------------------------------------------------------------------------------------
|	Fills in the `p' array with the pattern in column colN
*/
inline void NxsDiscreteMatrixBase::GetPattern(
  unsigned colN, 
  DataStorageType *p) const
  	{
  	const unsigned stateWidth = dataType.GetNumElementsPerCode();
  	for (unsigned rowN = 0; rowN < GetNumRows(); ++rowN)
  		{
  		const DataStorageType *toCopy = GetState(rowN, colN);
  		nxs_copy_n(toCopy, p, stateWidth);
  		p += stateWidth;
  		}
  	}
  	
/*--------------------------------------------------------------------------------------------------
|	Returns true if cell[rowN][colN] has all of its bit set to 1
*/
inline bool NxsDiscreteMatrixBase::IsMissing(
  unsigned rowN, 
  unsigned colN) const
  	{
  	return (ArrayCompare<DataStorageType>(GetState(rowN, colN), dataType.GetMissingBitCode(), nUnitsPerElement) == 0);
  	}

/*--------------------------------------------------------------------------------------------------
|	Returns true if cell[rowN][colN] has all of its bit set to 1
*/
inline bool NxsDiscreteMatrix::IsMissing(
  unsigned rowN, 
  unsigned colN) const
  	{
  	return IsGap(rowN, colN) && NxsDiscreteMatrixBase::IsMissing(rowN, colN);
  	}

/*--------------------------------------------------------------------------------------------------
|	returns a set of indices of columns in which the row specified by rowN has gaps
*/
inline const IndSet &  NxsDiscreteMatrix::GetGapMatrixRow(unsigned rowN) const
	{
	TwoDLookUp::const_iterator gIt = gapped.find(rowN);
  	if (gIt == gapped.end())
  		return emptyIndSet;
  	return gIt->second;
	}
	
/*--------------------------------------------------------------------------------------------------
|	returns true if cell[rowN][colN] was specified as a gap character
*/
inline bool NxsDiscreteMatrix::IsGap(
  unsigned rowN, 
  unsigned colN) const
  	{
  	assert(IsValidGetElement(rowN, colN));
	const IndSet & setOfGappedCells = GetGapMatrixRow(rowN);
  	return (setOfGappedCells.find(colN) != setOfGappedCells.end());
  	}
  	
/*--------------------------------------------------------------------------------------------------
|	returns true if cell[rowN][colN] was specified as a polymorphic character
*/
inline bool NxsDiscreteMatrixBase::IsPolymorphic(
  unsigned rowN, 
  unsigned colN) const
  	{
  	assert(IsValidGetElement(rowN, colN));
  	TwoDLookUp::const_iterator pIt = polymorphic.find(rowN);
  	if (pIt == polymorphic.end())
  		return false;
  	return (pIt->second.find(colN) != pIt->second.end());
  	}
  
/*--------------------------------------------------------------------------------------------------
|	Sets the cell[rowN][colN] as gap (zeroes the matrix element and adds the rowNxcolN to the map
|	of gapped sites)
*/
inline void NxsDiscreteMatrix::SetGap(
  unsigned rowN, 
  unsigned colN)
	{
	ZeroState(rowN, colN);
	AddGap(rowN, colN);
	}
	
/*--------------------------------------------------------------------------------------------------
|	Adds the rowNxcolN to the map of gapped sites and "zeros" the cell of the matrix
*/
inline void NxsDiscreteMatrix::AddGap(
  unsigned rowN, 
  unsigned colN)
	{
	gapped[rowN].insert(colN);
#	if defined(ZERO_WHEN_ADDING_GAP_STATES) 
		ZeroState(rowN, colN);
#	endif
	}
/*--------------------------------------------------------------------------------------------------
|	Sets the cell[rowN][colN] as missing (all 1's).
*/
inline void NxsDiscreteMatrixBase::SetMissing(
  unsigned rowN, 
  unsigned colN)
	{
	SetState(rowN, colN, dataType.GetMissingBitCode());
	}
/*--------------------------------------------------------------------------------------------------
|	Sets the cell[rowN][colN] as missing (all 1's).
*/
inline void NxsDiscreteMatrix::SetMissing(
  unsigned rowN, 
  unsigned colN)
	{
	AddGap(rowN, colN);
	SetState(rowN, colN, dataType.GetMissingBitCode());
	}
	
/*--------------------------------------------------------------------------------------------------
|	Adds the rowNxcolN to the map of polymorphic sites.
*/
inline void NxsDiscreteMatrixBase::SetPolymorphic(
  unsigned rowN, 
  unsigned colN)
	{
	polymorphic[rowN].insert(colN);
	}

/*--------------------------------------------------------------------------------------------------
|	returns true if the args describe a valid cell in the matrix whose state can be set.
*/
inline bool NxsUnalignedMatrix::IsValidPutElement(
  unsigned rowN,
  unsigned ) const
  	{
  	return (rowN < GetNumRows());
  	}
  	
/*--------------------------------------------------------------------------------------------------
|	adds indexOfState to cell[rowN][colN] 
*/
inline void NxsUnalignedMatrix::AddStateIndex(unsigned rowN, unsigned colN, unsigned indexOfState)
	{
	assert(IsValidPutElement(rowN, colN));
	PadRowZeroLast(rowN, colN);
	NxsDiscreteMatrixBase::AddStateIndex(rowN, colN, indexOfState);
	}
	
/*--------------------------------------------------------------------------------------------------
|	adds the bits in inValue to cell[rowN][colN]
*/
inline void NxsUnalignedMatrix::AddState(unsigned rowN, unsigned colN,const DataStorageType *inValue)
	{
	assert(IsValidPutElement(rowN, colN));
	PadRowZeroLast(rowN, colN);
	NxsDiscreteMatrixBase::AddState(rowN, colN, inValue);
	}
	
/*--------------------------------------------------------------------------------------------------
|	zeros cell[rowN][colN] 
*/
inline void NxsUnalignedMatrix::ZeroState(unsigned rowN, unsigned colN)
	{
	assert(IsValidPutElement(rowN, colN));
	if (PadRowZeroLast(rowN, colN))
		NxsDiscreteMatrixBase::ZeroState(rowN, colN);
	}

} // namespace ncl	
#endif
