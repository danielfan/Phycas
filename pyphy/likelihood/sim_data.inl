namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `divide_by' to 1.0, and `total_count' and `pattern_length' both to zero.
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
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a const reference to the data member `sim_pattern_map'.
*/
inline const PatternMapType & SimData::getSimPatternMap() const
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
	assert(ntaxa > 0);
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
|	Sets `tmp_pattern'[`pos'] to the supplied `state'.
*/
inline void SimData::setState(
  unsigned pos, 	/**< is the state to assign to every element in `tmp_pattern' */
  int8_t state)		/**< is the state to assign to every element in `tmp_pattern' */
	{
	assert(pos < pattern_length);
	tmp_pattern[pos] = state;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds `tmp_pattern' to `sim_pattern_map' then passes `missing_state' to the wipePattern() function to fill 
|	`tmp_pattern' with invalid values.
*/
inline void SimData::insertPattern(PatternCountType count)
	{
	// In debug build check to make sure there are no elements in tmp_pattern that still equal
	// missing_state (if assert trips, it means that not all elements of tmp_pattern have been 
	// replaced with actual states)
	assert(std::find(tmp_pattern.begin(), tmp_pattern.end(), missing_state) == tmp_pattern.end());

	// Add tmp_pattern to sim_pattern_map if it has not yet been seen, otherwise increment the count 
	// for this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
	PatternMapType::iterator lowb = sim_pattern_map.lower_bound(tmp_pattern);
	if (lowb != sim_pattern_map.end() && !(sim_pattern_map.key_comp()(tmp_pattern, lowb->first)))
		{
		// pattern is already in sim_pattern_map, increment count
		lowb->second += count;
		}
	else
		{
		// sim_pattern_map has not yet been stored in sim_pattern_map
		sim_pattern_map.insert(lowb, PatternMapType::value_type(tmp_pattern, count));
		}

	total_count += count;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds data currently stored in `sim_pattern_map' to the patterns already in `other'. The value of `mult' is used to
|	modify the counts before they are added to `other'; that is, the count of each pattern added to `other' is the 
|	original count multiplied by `mult'. Normally, `mult' would be specified to be 1.0, but in some cases it is
|	necessary to build up SimData objects that represent averages of other SimData objects, and it is these situations
|	where `mult' is handy. Assumes that `mult' is positive and non-zero. Assumes that `pattern_length' for this SimData
|	object is identical to the `pattern_length' of `other'.
*/
inline void SimData::addDataTo(SimData & other, PatternCountType mult)
	{
	assert(mult > 0.0);
	if (total_count == 0.0)
		return;
	if (other.getTotalCount() == 0)
		{
		other.resetPatternLength(pattern_length);
		}
	assert(pattern_length == other.getPatternLength());
	for (PatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
		{
		PatternCountType count = it->second;
		VecStateList & other_pattern = other.getCurrPattern();
		std::copy(it->first.begin(), it->first.end(), other_pattern.begin());
		other.insertPattern(mult*count);
		}
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
|	Returns the value of the `total_count' data member.
*/
inline PatternCountType SimData::getTotalCount()
	{
	return total_count;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Divides the count associated with every pattern in `sim_pattern_map' by `total'. Also divides `total_count' by
|	`total'. Assumes `total' is greater than zero. 
*/
inline void SimData::divideBy(PatternCountType total)
	{
	assert(total > 0.0);
	total_count /= total;
	for (PatternMapType::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
		{
		it->second /= total;
		}
	}

}	// namespace phycas
