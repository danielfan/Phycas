#ifndef NCL_NXS_STORED_MATRIX_H
#define NCL_NXS_STORED_MATRIX_H

#include "phypy/src/ncl/misc/nxs_discrete_matrix.hpp"
#include "phypy/src/ncl/characters/nxs_characters_manager.hpp"
#include "phypy/src/ncl/misc/eliminated_index_slider.hpp"

class NxsCharactersManager;
namespace ncl
{

class StoredMatrix
	{
	public:
		StoredMatrix(	unsigned nPrevChars, 
						DiscreteMatrixShPtr , 
						unsigned totalNChar, 
						const NxsIndexSet *, 
						const LabelMap *, 
						const StateLabelMap *);
		unsigned 		GetFirstGlobalIndex() const {return nCharsBefore;}
		unsigned 		size() const {return totalNChars;}
		
		void			GetPatternFromGlobalIndex(unsigned ind, DataStorageType *p) const
			{
			unsigned locInd = indexAdjuster.GetElementIndexFromGlobalIndex(ind);
			matrix->GetPattern(locInd, p);
			}
		DiscreteMatrixShPtr 	matrix;
		unsigned		getNCharsBefore() const
			{
			return nCharsBefore;
			}
		unsigned		getTotalNChars() const
			{
			return totalNChars;
			}
	private:
		void			IntializeNonMatrix(const NxsIndexSet *eliminated, const LabelMap *chLabels, const StateLabelMap *chStateLabels);
		
		unsigned		 		totalNChars;
		unsigned				nCharsBefore;
		LabelMap				charLabels;
		StateLabelMap			stateLabels;
		EliminatedIndexSlider 	indexAdjuster;
		
		friend class NxsCharactersManager;
	};

inline void StoredMatrix::IntializeNonMatrix(
  const NxsIndexSet *eliminated, 
  const LabelMap *chLabels, 
  const StateLabelMap *chStateLabels)
  	{
  	if (chLabels != NULL)
  		charLabels = *chLabels;
  	if (chStateLabels != NULL)
  		stateLabels = *chStateLabels;
  	indexAdjuster.SetNumPreviousIndices(nCharsBefore);
  	indexAdjuster.SetEliminated(*eliminated, totalNChars);
  	}
  	
/*--------------------------------------------------------------------------------------------------
|	Stores an alias to the matrix, and initializes the other fields.
*/
inline StoredMatrix::StoredMatrix(
  unsigned nChBefore,
  ncl::DiscreteMatrixShPtr mat, 
  unsigned totNChar, 
  const NxsIndexSet *eliminated,
  const LabelMap *chLabels, 
  const StateLabelMap *chStateLabels)
  	:matrix(mat),
  	totalNChars(totNChar),
  	nCharsBefore(nChBefore)
  	{
  	IntializeNonMatrix(eliminated, chLabels, chStateLabels);
  	}
} // namespace ncl
#endif
