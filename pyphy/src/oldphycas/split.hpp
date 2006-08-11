#ifndef PHO_SPLIT_H
#define PHO_SPLIT_H
#include <set>
#include <limits>
#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/misc/algorithm_extensions.hpp"
class TreeID;
class NxsTaxaManager;

typedef union
	{
	unsigned	i;
	double		f;
	} Score;

typedef Score			Length;
typedef unsigned long	split_t;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the notion of a taxon bipartition, or split. Each split is stored as a collection of bits, with the 
|	bits that are set representing the taxa above an edge in the tree. The type split_t defines the number of bits in 
|	one unit. If there are more taxa than bits in one unit, the split must be represented by multiple units. 
|	The STATIC_DATA_FUNC CalcStatics figures out the number of units that must be used.
*/
class Split
	{
	public:
							Split();
							Split(const NxsIndexSet & t);
							Split(const Split &other);
		
							~Split();
		
		STATIC_DATA_FUNC std::string 	GetStaticsInfo();
		STATIC_DATA_FUNC unsigned 		GetNTaxa() {return splitNTax;}
		STATIC_DATA_FUNC void 			CalcStatics(unsigned n);
		STATIC_DATA_FUNC double 	 	CalcLnPriorProb(unsigned n, unsigned m);
		
		unsigned			CalcComplexity() const;
		double				CalcLnPriorProb() const;
		double				CalcPriorProb() const;
		void	 			CreateAndAppendPatternRepresentation(std::string *) const;
		unsigned			CountOnBits() const;
		std::string			CreateIdRepresentation() const;
		std::string			CreateNewickRepresentation(const NxsTaxaManager *taxaMgr) const;
		std::string			CreatePatternRepresentation() const;
		void				GetIncludedExcludedRepresentation(NxsIndexSet *inc, NxsIndexSet *exc) const;
		bool 				Equals(const Split &other) const;
		bool 				IsBitSet(unsigned t) const;
		bool				IsCompatible(const Split &other) const;
		bool 				IsLessThan(const Split &other) const;
		bool 				SubsumedIn(const Split &other, unsigned startUnit = 0) const;

		bool				operator!=(const Split &other) const; //DANIEL 2-Dec-2004
		Split 				operator&(const Split &other) const; //DANIEL 2-Dec-2004
		Split 				operator|(const Split &other) const; //DANIEL 2-Dec-2004
		Split 				operator^(const Split &other) const; //DANIEL 2-Dec-2004
		Split 			  & operator^=(const Split &other); //DANIEL 2-Dec-2004
		Split			  & operator&=(const Split &other); //MTH 3-Dec-2004
		Split			  & operator|=(const Split &other); //MTH 3-Dec-2004

		bool				operator<(const Split &other) const;
		bool				operator==(const Split &other) const;
		Split 				&operator=(const Split &other);

		void 				Clear();
		void 				CombineWith(const Split &other);
		void 				IntersectWith(const Split &other);
		void 				SetBit(unsigned t);
		void 				SetBits(const NxsIndexSet & t);

		void 				UnsetBit(unsigned t); //DANIEL 2-Dec-2004

		void 				InvertSplit();
			
	private:
		void 				Resize();
		
	private:
		split_t					*unit;					/*< is the array of nunits split units */
		STATIC_CONST const unsigned		bits_per_unit = (CHAR_BIT) * sizeof(split_t);			/*< is the number of bits in a variable of type split_t */
		STATIC_CONST const split_t		unity = (split_t) 1 ;					/*< is a split_t variable with only the least significant bit set */
		STATIC_DATA unsigned		splitNTax;				/*< is the number of taxa currently under consideration */
		STATIC_DATA unsigned		nunits;					/*< is the length of the array necessary to represent a split */
		STATIC_DATA char			on_symbol;				/*< is the symbol used to represent bits that are set (i.e. "on") */
		STATIC_DATA char			off_symbol;				/*< is the symbol used to represent bits that have been cleared (i.e. "off") */
		STATIC_DATA split_t			mask;					/*< a split_t variable used to exclude the remainder bits at the end of the final unit when `splitNTax' modulo `bits_per_unit' is not 0 */

		friend std::istream &operator>>(std::istream &in, TreeID &id);
		friend std::istream &operator>>(std::istream &in, Split &s);
		friend std::string &operator<<( std::string &in, const TreeID &id);
		friend std::string &operator<<( std::string &in, const Split &s);
		friend class Tree;
		friend class SplitManager;
	};

typedef std::set<Split> SplitSet;

unsigned FindSplitsAbsentInTestTree(const SplitSet &refTree, const SplitSet &testTree, SplitSet *missing);

//inline STATIC_DATA unsigned Split::GetNTax()
//	{
//	return ntax;
//	}

// Daniel 2-Dec-2004
/*----------------------------------------------------------------------------------------------------------------------*/
inline bool Split::operator!=(const Split &other) const
	{
		if (*this == other)
			return false;
		else
            return true;
	}
	

// Daniel 2-Dec-2004 MTH converted operator& to a call to operator&=
/*----------------------------------------------------------------------------------------------------------------------*/
inline Split Split::operator&(const Split &other) const
	{
	Split tmp(*this);
	return tmp &= other;
	}
	
/*----------------------------------------------------------------------------------------------------------------------*/
inline Split & Split::operator&=(const Split &other)
	{
	for (unsigned i = 0; i < nunits; ++i)
		unit[i] &= other.unit[i];
	return *this;
	}

// Daniel 3-March-2005
inline Split Split::operator^(const Split &other) const
        {
        Split tmp(*this);
        return tmp ^= other;
        }

// Daniel 3-March-2005
inline Split & Split::operator^=(const Split &other)
        {
        for (unsigned i = 0; i < nunits; ++i)
                unit[i] ^= other.unit[i];
        return *this;
        }

// Daniel 2-Dec-2004 MTH converted operator! to a call to operator|=
/*----------------------------------------------------------------------------------------------------------------------*/
inline Split Split::operator|(const Split &other) const
	{
	Split tmp(*this);
	return tmp |= other;
	}

/*----------------------------------------------------------------------------------------------------------------------*/
inline Split & Split::operator|=(const Split &other)
	{
	for (unsigned i = 0; i < nunits; ++i)
		unit[i] |= other.unit[i];
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls IsLessThan(`other') and returns the value returned by that function.
*/
inline bool Split::operator<(const Split &other) const
	{
	return IsLessThan(other);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Equals(`other') and returns the value returned by that function.
*/
inline bool Split::operator==(const Split &other) const
	{
	return Equals(other);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates unit array nunits long, copies values from other's unit array, then returns reference to *this.
*/
inline Split &Split::operator=(const Split &other)
	{
	nxs_copy(other.unit, other.unit + nunits, unit);
	return *this;
	}



/*----------------------------------------------------------------------------------------------------------------------
|	Creates unit array nunits long and then copies values from other's unit array using operator=.
*/
inline Split::Split(const Split& other)
	{
	unit = new split_t[nunits];
	*this = other;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes the unit array.
*/
inline Split::~Split()
	{
	delete [] unit;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates unit array nunits long and calls Clear().
*/
inline Split::Split()
	{
	unit = new split_t[nunits];
	Clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the prior probability of this split given flat prior on topologies. Specifically, returns exp(x), where x 
|	is the result from CalcLnPriorProb().
*/
inline double Split::CalcPriorProb() const
	{
	return std::exp(CalcLnPriorProb()); //@ why is the std:: necessary?
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the log of the prior probability of this split given flat prior on topologies. Specifically, returns the 
|	result of calling CalcLnPriorProb(n, m), where n is the number of set bits and m the number of cleared bits, with
|	n + m equal to the number of taxa.
*/
inline double Split::CalcLnPriorProb() const
	{
	// n is the number of bits set
	// m is the number not set
	// n + m should add up to Split::ntax
	//
	const unsigned n = CountOnBits();
	const unsigned m = Split::splitNTax - n;
	return CalcLnPriorProb(n, m);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if all elements of the unit array are equal to all elements of the Split object other.
*/
inline bool Split::Equals(const Split& other) const
	{
	//@POL what about irrelevant bits (i.e. if only 5 taxa, then 3 bits must go unused, are these guaranteed to 
	// all be set the same way in different Splits?)
	return std::equal(unit, unit + nunits, other.unit);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets all nunits values in the unit array to 0.
*/
inline void Split::Clear()
	{
	replace_all(unit, unit + nunits, (split_t)0);
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Performs a bitwise OR of each element of unit with the corresponding element in other's unit array. The result is
|	a union of the sets defined by the two Split objects, and is useful in creating the Split for an interior node, which
|	is the union of the splits of its immediate descendants.
*/
inline void Split::CombineWith(const Split &other)
	{
	*this |= other;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts this split into the intersection between this split and the other split. The other split remains unaffected
|	but this split object is modified such that each element in the unit array undergoes a bitwise AND operation with
|	the corresponding element from other's unit array.
*/
inline void Split::IntersectWith(
	const Split & other)	/* the other split */
	{
	*this &= other;
	}

#endif
