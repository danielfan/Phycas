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
#include <fstream>
namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `total_count' and `pattern_length' both to zero.
*/
inline SimData::SimData()
  : total_count(0.0), pattern_length(0)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets this object to its just-constructed state.
*/
inline void SimData::clear()
	{
	tmp_pattern.clear();
	sim_pattern_map.clear();
	total_count = 0.0;
	pattern_length = 0;
#if POLPY_NEWWAY
    patternVect.clear();
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a const reference to the data member `sim_pattern_map'.
*/
inline const SimPatternMapType & SimData::getSimPatternMap() const
	{
	return sim_pattern_map;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `ntaxa' equals `pattern_length', this function returns immediately without taking any action. If `ntaxa' differs
|	from `pattern_length', sets `pattern_length' to `ntaxa' and clears `sim_pattern_map' because it is invalidated when 
|	`pattern_length' changes. Also refreshes the `tmp_pattern' data member to conform to the new value of 
|	`pattern_length'. Assumes `ntaxa' is non-zero. Use clear() to set `ntaxa' to 0 and return the object to its 
|	just-constructed state.
*/
inline void SimData::resetPatternLength(
  unsigned ntaxa)	/**< is the number of taxa (same as the number of elements in a pattern vector) */
	{
	PHYCAS_ASSERT(ntaxa > 0);
	if (ntaxa == pattern_length)
		return;

	clear();
	pattern_length = ntaxa;

	// Create a VecStateList vector with ntaxa elements all of which are -1
	// and swap into tmp_pattern so that tmp_pattern now has the correct size
	VecStateList v(ntaxa, missing_state);	// GCC 3.2.3 (Red Hat Linux 3.2.3-20) requires this split 
	tmp_pattern.swap(v);					// into two lines 
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets every element in `tmp_pattern' to the value of `missing_state'.
*/
inline void SimData::wipePattern()
	{
	tmp_pattern.assign((VecStateList::size_type)pattern_length, SimData::missing_state); 
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets all counts in `sim_pattern_map' to zero.
*/
inline void SimData::zeroCounts()
	{
	for (SimPatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
		{
		it->second = 0.0;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Saves all counts as a tab-delimited row appended to the file whose name is 'fn'. Note: no header is output
|	identifying which pattern goes with each count, so this will not be useful unless only the unidentified counts are 
|	of interest.
*/
inline void SimData::appendCountsToFile(std::string fn, bool binary)
	{
	std::ofstream outf(fn.c_str(), std::ios::out | std::ios::app);
	SimPatternMapType::iterator it = sim_pattern_map.begin();
	unsigned count = (unsigned)(it->second);
	if (binary && count > 0)
		count = 1;
	outf << count;
	++it;
	for (; it != sim_pattern_map.end(); ++it)
		{
		count = (unsigned)(it->second);
		if (binary && count > 0)
			count = 1;
		outf << '\t' << count;
		}
	outf << std::endl;
	outf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Saves all counts as a tab-delimited row appended to the file whose name is 'fn'. Note: no header is output
|	identifying which pattern goes with each count, so this will not be useful unless only the unidentified counts are 
|	of interest.
*/
inline void SimData::debugAppendCountsToFile(std::string row_name, std::string fn)
	{
	std::ofstream outf(fn.c_str(), std::ios::out | std::ios::app);
	outf << row_name;

#if 0

	for (unsigned i = 0; i < 4; ++i)
		{
		for (unsigned j = 0; j < 4; ++j)
			{
			for (unsigned k = 0; k < 4; ++k)
				{
					for (unsigned m = 0; m < 4; ++m)
						{
						VecStateList v;
						v.push_back(i);
						v.push_back(j);
						v.push_back(k);
						v.push_back(m);
						SimPatternMapType::iterator it = sim_pattern_map.find(v);
						PatternCountType count = 0.0;
						if (it != sim_pattern_map.end())
							{
							count = (PatternCountType)(it->second);
							}
						outf << '\t' << count;
						}
				}
			}
		}

#else

	SimPatternMapType::iterator it = sim_pattern_map.begin();
	PatternCountType count = (PatternCountType)(it->second);
	outf << count;
	++it;
	for (; it != sim_pattern_map.end(); ++it)
		{
		count = (PatternCountType)(it->second);
		outf << '\t' << count;
		}

#endif

	outf << std::endl;
	outf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of patterns, each represented as a string. The patterns are in the same order as the counts 
|	returned by the appendCountsToFile function, so this function can be used to identify the patterns belonging to the
|	pattern counts returned by that function. This is a slow function, so don't use it in situations where speed is 
|	critical.
*/
inline std::vector<std::string> SimData::getPatterns(std::vector<std::string> symbols)
	{
	std::vector<std::string> v;
	for (SimPatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
		{
		// Create a string out of the pattern
		std::string s;
		for (VecStateList::const_iterator i = it->first.begin(); i != it->first.end(); ++i)
			{
			s += symbols.at(*i);
			}
		v.push_back(s);
		}
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Saves all counts as a string that can be used within Maple to create a 3D surface plot. To make a surface plot of 
|	the following 2x3 matrix in Maple 9,
|>
|	1  2  3
|	4  6  8
|>
|	place the following commands in a file named maple_commands:
|>
|	with(linalg);
|	with(plots);
|	points := [[[0,0,1],[0,1,2],[0,2,3]],[[1,0,4],[1,1,6],[1,2,8]]];
|	surfdata(points, axes=framed, labels=["rows", "cols", "height"]);
|>
|	then read this file by typing 'read maple_commands;' in the Maple interpreter. In our case, rows will be separate
|	posterior predictive data sets and columns will be patterns (256 columns for 4 taxa and DNA data). The height
|	will be the number of counts. This function returns a string representing one row of the points matrix. For 
|	example, if this data set were going to be the 2nd row (index 1) in the above example, the string returned would be 
|>
|	[[1,0,4],[1,1,5],[1,2,6]]
|>
|	
*/
inline std::string SimData::createMapleTuples(unsigned row, unsigned cutoff)
	{
	bool use_cutoff = true;
	if (cutoff == 0)
		use_cutoff = false;

	std::string s = "[";
	unsigned col = 0;
	SimPatternMapType::iterator it = sim_pattern_map.begin();
	unsigned count = (unsigned)(it->second);
	if (use_cutoff && count > cutoff)
		count = 0;
	double log_count = count > 0 ? std::log((double)count) : 0.0;
	s += str(boost::format("[%d,%d,%f]") % row % col % log_count);
	++it;
	++col;
	for (; it != sim_pattern_map.end(); ++it, ++col)
		{
		count = (unsigned)(it->second);
		if (use_cutoff && count > cutoff)
			count = 0;
		log_count = count > 0 ? std::log((double)count) : 0.0;
		s += str(boost::format(",[%d,%d,%f]") % row % col % log_count);
		}
	s += "]";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `tmp_pattern'[`pos'] to the supplied `state'.
*/
inline void SimData::setState(
  unsigned pos, 	/**< is the position in `tmp_pattern' to set */
  int8_t state)		/**< is the state to assign to the element in `tmp_pattern' at position pos */
	{
	PHYCAS_ASSERT(pos < pattern_length);
	tmp_pattern[pos] = state;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the `tmp_pattern' data member.
*/
inline VecStateList & SimData::getCurrPattern()
	{
	return tmp_pattern;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the `pattern_length' data member.
*/
inline unsigned SimData::getPatternLength()
	{
	return pattern_length;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the `total_count' data member, which represents the total number of patterns added (not 
|	`sim_pattern_map'.size()).
*/
inline PatternCountType SimData::getTotalCount()
	{
	return total_count;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of the `total_count' data member, which represents the sum of all pattern counts (not 
|	`sim_pattern_map'.size()).
*/
inline void SimData::setTotalCount(PatternCountType total)
	{
	total_count = total;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value `sim_pattern_map'.size(), the number of unique patterns currently stored.
*/
inline unsigned SimData::getNUniquePatterns()
	{
	return (unsigned)sim_pattern_map.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Divides the count associated with every pattern in `sim_pattern_map' by `factor'. Also divides `total_count' by
|	`factor'. Assumes `factor' is greater than zero. 
*/
inline void SimData::divideBy(PatternCountType factor)
	{
	PHYCAS_ASSERT(factor > 0.0);

	total_count /= factor;

	for (SimPatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
		{
		it->second /= factor;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Multiplies the count associated with every pattern in `sim_pattern_map' by `factor'. Also multiplies `total_count'
|	`factor'. Assumes `factor' is greater than zero. 
*/
inline void SimData::multBy(PatternCountType factor)
	{
	if (sim_pattern_map.empty())
		return;

	total_count *= factor;

	for (SimPatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
		{
		it->second *= factor;
		}
	}
}	// namespace phycas
