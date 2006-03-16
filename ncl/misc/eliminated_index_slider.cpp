#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/misc/eliminated_index_slider.hpp"


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
		assert(prevElement <= maxIndex);
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
				assert(ind - currIndex == mIt->second.second);
			for (; ind <= mIt->first;)	
				{
				v->push_back(UINT_MAX);
				++ind;			//also iterating ind here 
				if (ind >= maxSize)
					return;
				}
			++mIt;
			assert(mIt != offsets.end());
			begOfNextBlock = (mIt->first - mIt->second.first);
			assert(ind < begOfNextBlock);

			}
		v->push_back(currIndex++);
		}
	}

