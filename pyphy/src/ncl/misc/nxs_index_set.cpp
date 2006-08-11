//#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/misc/nxs_index_set.hpp"
#include "ncl/misc/string_extensions.hpp"
using std::string;

/*-------------------------------------------------------------------------------------------------------------------------- 
| 	Appends the nexus description of the set s (e.g "1-10 \2 15 21" )
*/
string & operator+= (string & st, const NxsIndexSet &iSet)
	{
	st += iSet.GetNexusDescription();
	return st;
	}
	
string NxsIndexSet::GetNexusDescription(
  bool addOne) const	/*add one for output to user */
	{
	string retStr;
	if (size() == 0)
		return retStr;
	
	const_iterator c = begin();
	unsigned addVal = (addOne ? 1U : 0U);
	unsigned begRange = addVal + *c++;
	retStr << begRange;
	if (c == end()) // only one index in the set
		return retStr;
	retStr << ' ';
	unsigned prevC = addVal + *c++;
	unsigned currC = prevC;
	unsigned stride = prevC - begRange;
	if (c == end())
		{
		retStr << prevC;
		return retStr;
		}
	while (c != end())
		{
		currC = addVal + *c++;
		unsigned nStride = currC - prevC;
		if (nStride == stride)	
			{
			// continuing the range e.g. beg = 1 prev = 3 curr = 5 then we just want to set prev to 5
			prevC = currC;
			}
		else
			{
			if (prevC - begRange == stride)
				{	
				//	not really a range with an even stride, just three values 
				//	e.g. beg = 1 prev = 3 curr = 10 then we want to 
				//	set beg to 3 and add it to the string and set the stride to 7 and prev to 10
				//
				retStr << prevC << ' ';
				begRange = prevC;
				prevC = currC;
				stride = nStride;
				}
			else
				{
				//	prev is more than on stride from the beginning spit out a range
				//	e.g. beg = 1 [then 3] prev = 5 curr = 10 then we want to 
				//	add "-5 \2 10" to the string set beg to 10 and read in a new value (if possible for prev)
				retStr << "- " << prevC;
				if (stride > 1)
					retStr << '\\' << stride;
				begRange = prevC = currC;
				retStr << ' ' << begRange << ' ';
				if (c != end())
					{
					prevC = addVal + *c++;
					stride = prevC - begRange;
					}
				}
			}
		}
	if (currC != begRange)
		{
		//	we still need to spit out the value in prevC
		//
		if (prevC - begRange == stride)
			retStr << prevC;
		else
			{
			retStr << "- " << prevC;
			if (stride > 1)
				retStr << '\\' << stride << ' ';
			}
		}
	return retStr;
	}

void NxsIndexSet::Mask(
  const NxsIndexSet &m)
  	{
  	NxsIndexSet temp;
  	temp.SetToIntersection(*this, m);
  	*this = temp;
  	}
  	
/*-------------------------------------------------------------------------------------------------------------------------- 
| 	Walks passed n elements and returns the next one (or UINT_MAX is n >= size())
*/
unsigned NxsIndexSet::FindNthElement(unsigned n) const
	{
	NxsIndexSet::const_iterator nIt = begin();
	for (unsigned i = 0; nIt != end() && i < n; ++i)
		++nIt;
	if (nIt == end())
		return UINT_MAX;
	return *nIt;
	}

/*-------------------------------------------------------------------------------------------------------------------------- 
| 	Counts the number of elements before matching is reached (returns UINT_MAX is matching isn't found)
*/
unsigned  NxsIndexSet::CountIndicesBeforeMember(unsigned matching) const
	{
	
	NxsIndexSet::const_iterator nIt = begin();
	for (unsigned i = 0; nIt != end(); ++i, ++nIt)
		{
		if (*nIt == matching)
			return i;
		}
	return UINT_MAX;
	}

void NxsIndexSet::AddToAllIndices(
  int n)
	{
	if (n == 0)
		return;
	const NxsIndexSet temp = *this;
	clear();
	NxsIndexSet::const_iterator temIt = temp.begin();
	if (n < 0)
		{
		n = -n;
		unsigned nPos = (unsigned) n;
		for (; temIt != temp.begin(); ++temIt)
			{
			if (*temIt >= nPos)
				insert(*temIt - nPos);
			}
		}
	else
		{
		for (; temIt != temp.begin(); ++temIt)
			insert(*temIt + n);
		}
	}

void NxsIndexSet::Crop(
  unsigned lowest, 
  unsigned highest)
  	{
  	assert(highest >= lowest);
  	const unsigned currLast = GetLast();
  	const unsigned currFirst = GetFirst();
  	if (lowest > currLast || highest < currFirst)
  		{
  		clear();
  		return;
  		}
  	if (lowest > 0 && lowest - 1 >= currFirst)
		EraseRange(currFirst, lowest - 1);
	if (highest + 1 <= currLast)
		EraseRange(highest + 1, currLast);
	}	

