#ifndef COMPRESSED_2D_ARRAY_H
#define COMPRESSED_2D_ARRAY_H

/*----------------------------------------------------------------------------------------------------------------------
| 	This structure is for creating cache-friendly 2D arrays when the length of the arrays varies.  
|	This class holds arrays of type nStateInt and was designed for storing data at the tips of trees.  In these cases 
|	there is usually one state seen (0, 1, 2, or 3) but occasionally there will be more.  So if we want to 
|	store the array of states :
|	0 3 1 1 {2,3} 0	{0,1,2,3} 2
|	we could use a 2-D array States[i][j] where i indexes the site and j indexes the state.  To cover all bases we'd 
|	have to make the length of the array 4 * TheNumberOfSites.  This is wasteful and will slow down likelihood 
|	calculations. Instead we store them in a Compressed2DArray as the array:
|	0 3 1 1 -2 2 3 0 -4 2
|	To get the correct sites out of this array you have to walk down the array following the rules:
|	if ( *arr >= 0)  then *arr is the value (what would be indexed as States[i][0] in the traditional scheme)
|	if ( *arr < 0)	then - *arr is the number of states in this "subarray"  the states will then follow.
|	if ( *arr == - maximumNumber)	this indicates "complete ambiguity" - all of the states are present they will NOT
|	be listed afterwards (since you know what they will be).
|
|	This creates a very compact, storage that is easy to interpret and linear (only one array has to be accessed).
|	The main disadvantage is that if you need the states at site 53, you have to walk through the array (since how far
|	that site will be from the beginning depends on the amount of ambiguity).
|
|	This notational style is used for other frequently-accessed arrays of arrays of positive, bounded integers.
*/
template <typename T> class Compressed2DArray
	{
	public :
		unsigned nStates;//POL-121803 
		std::vector<T> arr;

		//@POL Mark, I am a little worried about assigning a default value of 0 for nStates in the constructor
		// Need to talk about where you use Compressed2DArray - can you use OrdCodedArrs in those places, given
		// that the only functions of Compressed2DArray (Append and AppendEmpty) require knowledge of the number
		// of states
		//

		//POL-121803 Compressed2DArray(unsigned lenGuess = 0) 
		Compressed2DArray(unsigned lenGuess = 0, unsigned numStates = 0) : nStates(numStates)
			{
			// Compressed2DArray should only be used with signed types. Would use a compile-time assert here 
			// if numeric_limits was more portable
			//
			PHYCAS_ASSERT(0 > (T) (-(T)1)); 

			if (lenGuess != 0)
				arr.reserve(lenGuess);
			}

		//POL-121803 Compressed2DArray(const vector<VecUInt> &full2DVec);
		Compressed2DArray(const VecVecUInt & full2DVec, unsigned numStates);

		void	Append(const VecUInt &);
		//POL-121803 void 	AppendEmpty(unsigned nStates)
		void 	AppendEmpty()
			{
			NXS_ASSERT(nStates != 0); // should not be using AppendEmpty if object instantiated with default constructor
			arr.push_back((T)(-(T)nStates));
			}
	};
	

/*----------------------------------------------------------------------------------------------------------------------
| 	Appends vector that represents one "row" in the flattened matrix appends to the flattened array vec
*/
template <typename T>
 void Compressed2DArray<T>::Append(
	const VecUInt &v) 
  	{
  	//NXS_ASSERT(!v.empty()); //@POL Mark, until we find a better solution, treating gaps as missing
	NXS_ASSERT(v.size() < 8*sizeof(T));
	unsigned sz = (unsigned)v.size();
	if (sz == 0)
		{
		// Gaps cause v to be empty
		//
		arr.push_back((T)(-(T)nStates));
		}
	else if (sz == 1)
		{
		// v contains just one state (not a gap, not missing, and no ambiguity or polymorphism)
		//
		arr.push_back((T) v[0]);
		}
	else
		{
		arr.push_back((T) (-(T)sz));
		if (sz < nStates) //POL-121803
			{
			for (VecUInt::const_iterator elIt = v.begin(); elIt != v.end(); ++elIt)
				arr.push_back((T) *elIt);
			}
		}
  	}	
 
/*----------------------------------------------------------------------------------------------------------------------
| 	
*/
template <typename T>
inline Compressed2DArray<T>::Compressed2DArray(
	const VecVecUInt &full2DVec,
	unsigned numStates) 
	: nStates(numStates) 
  	{
  	NXS_ASSERT(0 > (T) (-(T)1)); // Compressed2DArray only be used with signed types.  would use a compile-time assert if numeric_limits was more portable
	arr.reserve(full2DVec.size());
  	for (VecVecUInt_ConstIt fIt = full2DVec.begin(); fIt != full2DVec.end(); ++fIt)
  		Append(*fIt);
  	}

/*----------------------------------------------------------------------------------------------------------------------
|	OrdCodedArrs is derived from Compressed2DArray, differing mainly in having a specialized 
|	const_iterator that can handle ambiguity in state for a particular taxon. Here is how 
|	characters are stored in an Compressed2DArray object:
|>
|	States for 8 taxa for one pattern: 0 3 1 1 {2,3} 0 {0,1,2,3} 2
|	  (note ambiguity for taxa 5 and 7)
|	
|	As stored in Compressed2DArray:    0 3 1 1 -2 2 3 0 -4 2
|	  -2 indicates ambiguity/polymorphism and next 2 states are for same taxon
|	  -4 indicates completely missing data and there is no need to list states afterwards
|<
|	The OrdCodedArrs::const_iterator would return the following vectors for each of the 
|	eight taxa if the array is walked using DerefThenAdvance:
|<
|	(0)    0 is the only state for taxon 0
|	(3)    3 is the only state for taxon 1
|	(1)    1 is the only state for taxon 2
|	(1)    1 is the only state for taxon 3
|	(2,3)  taxon 4 is ambiguous, either state 2 or state 3 is present
|	(0)    0 is the only state for taxon 5
|	()     taxon 6 has completely missing data, so empty vector returned
|	(2)    2 is the only state for taxon 7
|<
*/
class OrdCodedArrs : public Compressed2DArray<NStateInt>
	{
	public :
		class const_iterator 
			{
			typedef std::vector<NStateInt>::const_iterator SubIterator;
			SubIterator vecIterator;
			unsigned nStates;

			public:
				bool operator==(const const_iterator r) const
					{
					return (vecIterator == r.vecIterator);
					}

				bool operator!=(const const_iterator r)	const
					{
					return !((*this) == r);
					}

				/** Equivalent to *iterator++ in a std::vector, except here dereferencing requires translating the subarray if there is ambiguity. */
				std::vector<unsigned> DerefThenAdvance()
					{
					std::vector<unsigned> retVec;
					const NStateInt el = *vecIterator++;
					if (el >= 0)
						retVec.push_back((unsigned) el);
					else
						{
						unsigned uEl = (unsigned) (-el);
						if (uEl != nStates)
							{
							for (unsigned i = 0; i < uEl; ++i, ++vecIterator)
								retVec.push_back((unsigned) *vecIterator);
							}
						}
					return retVec;
					}	

				/** Returns a vector of the states in the OrdCodedArrs element. Returns an empty vector for the all-states-ambiguous code. */
				std::vector<unsigned> operator*() const
					{
					std::vector<unsigned> retVec;
					const NStateInt el = *vecIterator;
					if (el >= 0)
						retVec.push_back((unsigned) el);
					else
						{
						unsigned uEl = (unsigned) (-el);
						// note that for completely missing data, the returned retVec is empty
						if (uEl != nStates)
							{
							SubIterator temp = vecIterator;
							++temp;
							for (unsigned i = 0; i < uEl; ++i, ++temp)
								retVec.push_back((unsigned) *temp);
							}
						}
					return retVec;
					}	
					
				const_iterator operator++() //prefix
					{
					const NStateInt el = *vecIterator;
					if (el >= 0)
						++vecIterator;
					else
						{
						unsigned uEl = (unsigned) (-el);
						if (uEl == nStates)
							++vecIterator;
						else
							std::advance<SubIterator>(vecIterator, 1 + uEl);
						}
					return *this;
					}

				const_iterator operator++(int)//postfix
					{
					const_iterator ret(*this); 
					++(*this);
					return ret;
					}	

				const_iterator(const SubIterator &r, unsigned nS)
				  : vecIterator(r), nStates(nS)
					{
					}					
			};
			
		OrdCodedArrs(unsigned numStates) : Compressed2DArray<NStateInt>(0, numStates) {}
		OrdCodedArrs(const VecVecUInt & v, unsigned numStates): Compressed2DArray<NStateInt>(v, numStates) {}
		OrdCodedArrs(const std::vector<DataStorageType *> &, unsigned wordLen, unsigned numStates); // see "constsiteinfo.cpp" for definition
		
		const_iterator begin() const
			{
			return const_iterator(arr.begin(), nStates);
			}
	};
	
//POL-121803 removed unsigned nStates; data member
//POL-121803 removed OrdCodedArrs(unsigned numStates): Compressed2DArray<NStateInt>(), nStates(numStates) {}
//POL-121803 removed OrdCodedArrs(const vector<VecUInt> &v, unsigned numStates): Compressed2DArray<NStateInt>(v), nStates(numStates) {}
//POL-121803 removed OrdCodedArrs(const vector<DataStorageType *> &, unsigned wordLen, unsigned numStates);
		
#endif
