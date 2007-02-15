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

#ifndef NCL_NXS_STORED_MATRIX_H
#define NCL_NXS_STORED_MATRIX_H

#include "phycas/src/ncl/misc/nxs_discrete_matrix.hpp"
#include "phycas/src/ncl/characters/nxs_characters_manager.hpp"
#include "phycas/src/ncl/misc/eliminated_index_slider.hpp"

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
