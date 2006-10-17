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
#include "phypy/src/ncl/misc/eliminated_index_slider.hpp"


/*--------------------------------------------------------------------------------------------------
|	
*/
void EliminatedIndexSlider::SetEliminated(const NxsIndexSet &s, unsigned maxIndex)
	{
	typedef std::pair<unsigned, unsigned> PairUInt;
	clear();
	if (s.empty())
		return;
	
	NxsIndexSet::const_iterator setIt = s.begin();
	unsigned beginingOfBlock = *setIt;
	unsigned prevElement = beginingOfBlock;
	unsigned totalOffset = 0;
	
	for (setIt++ ; setIt != s.end(); ++setIt)
		{
		if (*setIt != prevElement + 1)
			{
			const unsigned sizeOfBlockMinusOne = beginingOfBlock- prevElement; 
			offsets[beginingOfBlock + sizeOfBlockMinusOne] = PairUInt(sizeOfBlockMinusOne, totalOffset);
			totalOffset += sizeOfBlockMinusOne + 1;
			prevElement = beginingOfBlock = *setIt;
			}
		else
			++prevElement;
		NXS_ASSERT(prevElement <= maxIndex);
		}
	const unsigned secSizeOfBlockMinusOne = beginingOfBlock- prevElement; 
	offsets[beginingOfBlock + secSizeOfBlockMinusOne] = PairUInt(secSizeOfBlockMinusOne, totalOffset);
	totalOffset += secSizeOfBlockMinusOne + 1;
	
	//	now add a block of length 1 at maxIndex + 1 which will provide the offset for the tail of the indices
	//
	offsets[maxIndex + 1] = PairUInt(0, totalOffset);
	}
	
void EliminatedIndexSlider::BuildOldToNewIndexMap(std::vector<unsigned> *v, unsigned maxSize) const
	{
	v->reserve(v->size() + maxSize);
	if (offsets.empty())
		{
		for (unsigned i = 0; i < maxSize; ++i)
			v->push_back(i);
		return;
		}
	OffsetMap::const_iterator mIt = offsets.begin();
	unsigned begOfNextBlock = (mIt->first - mIt->second.first);
	unsigned currIndex = 0;
	for (unsigned ind = 0; ind < maxSize; ++ind)
		{
		if (ind  >= begOfNextBlock)
			{
				NXS_ASSERT(ind - currIndex == mIt->second.second);
			for (; ind <= mIt->first;)	
				{
				v->push_back(UINT_MAX);
				++ind;			//also iterating ind here 
				if (ind >= maxSize)
					return;
				}
			++mIt;
			NXS_ASSERT(mIt != offsets.end());
			begOfNextBlock = (mIt->first - mIt->second.first);
			NXS_ASSERT(ind < begOfNextBlock);

			}
		v->push_back(currIndex++);
		}
	}

