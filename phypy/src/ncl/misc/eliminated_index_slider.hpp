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

#ifndef NCL_ELIMINATED_INDEX_SLIDER_H
#define NCL_ELIMINATED_INDEX_SLIDER_H

#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/misc/nxs_index_set.hpp"

/*--------------------------------------------------------------------------------------------------
|	A class that translates one indexing system into another.  
|	Designed for helping find characters after a characters block was read.  
|	When a new characters block is read its indeces are "off" by PREV_INDICES (the characters block 
|	starts numbering at 0, but this is really character PREV_INDICES in the global numbering)
|	The block has an index set of characters that were skipped
*/
class EliminatedIndexSlider
	{
	public:
		EliminatedIndexSlider() : nPreviousIndices(0){}
		void		clear();
		unsigned	GetElementIndexFromLocalIndex(unsigned) const;
		unsigned	GetElementIndexFromGlobalIndex(unsigned) const;
		void		SetNumPreviousIndices(unsigned n) { nPreviousIndices = n;}
		void		SetEliminated(const NxsIndexSet &, unsigned maxIndex);
		void		BuildOldToNewIndexMap(std::vector<unsigned> *v, unsigned) const;
		
	private:
		unsigned	UseOffsetMapToSlideIndex(unsigned) const;
		
		
		typedef std::pair<unsigned, unsigned> OffsetInfo;	//size of current elimated block - 1, offset BEFORE the block
		typedef std::map<unsigned, OffsetInfo> OffsetMap;
		OffsetMap 	offsets;
		unsigned	nPreviousIndices;
	};


inline unsigned EliminatedIndexSlider::GetElementIndexFromLocalIndex(unsigned localInd) const
	{
	if (offsets.empty()) 
		return localInd;
	return UseOffsetMapToSlideIndex(localInd);
	}

inline unsigned EliminatedIndexSlider::GetElementIndexFromGlobalIndex(unsigned globInd) const
	{
	globInd -= nPreviousIndices;
	return GetElementIndexFromLocalIndex(globInd);
	}

inline void EliminatedIndexSlider::clear()
	{
	offsets.clear();
	}

inline unsigned EliminatedIndexSlider::UseOffsetMapToSlideIndex(unsigned orig) const
	{
	OffsetMap::const_iterator mIt = offsets.lower_bound(orig);
	PHYCAS_ASSERT(mIt != offsets.end());  // as long as orig <= maxIndex sent in AddEliminated we should find an element.
#	if 0
		const unsigned &endOfEliminatedBlock = mIt->first;	
		const unsigned &sizeOfEliminatedBlock =  mIt->second.first;
		if (orig > endOfEliminatedBlock - sizeOfEliminatedBlock)
			return UINT_MAX;	//eliminated
		const unsigned &offsetBeforeEliminatedBlock = mIt->second.second;	
		return orig - offsetBeforeEliminatedBlock;
#	else	//same code, but less readable

		if (orig  >= (mIt->first - mIt->second.first))
			return UINT_MAX;	//eliminated
		return orig - mIt->second.second;

#	endif
	}

	
#endif

