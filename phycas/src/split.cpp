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

#include <cmath>
#include <set>
#include <limits>
#include <algorithm>
#include <boost/format.hpp>
#include "phycas/src/split.hpp"
#include "phycas/src/xphylogeny.hpp"

#define SPLIT_UNITY_VALUE   ((split_t)1)
#define BITS_PER_UNIT_VALUE ((CHAR_BIT)*sizeof(split_t))

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Clear().
*/
Split::Split()
  : bits_per_unit(BITS_PER_UNIT_VALUE), split_unity(SPLIT_UNITY_VALUE)
	{
	Clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `split_ntax' to 4, `nunits' to 1, `on_symbol' to an asterisk, `off_symbol' to hyphen, `excl_symbol' to x, 
|   `mask' to 0 and resizes the `units' vector to have a length of `nunits' and fills it with 0s.
*/
void Split::Clear()
	{
    split_ntax		= 4;
    nunits			= 1;
    on_symbol		= '*';
    off_symbol		= '-';
    excl_symbol		= 'x';
    mask			= (split_t)0;
	unit.resize(nunits, (split_t)0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a copy of the supplied `other' Split object by calling operator=.
*/
Split::Split(
  const Split & other) /**< is the Split object to be copied */
  : bits_per_unit(BITS_PER_UNIT_VALUE), split_unity(SPLIT_UNITY_VALUE)
	{
	*this = other;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a copy of the Split object `other', then returns reference to *this.
*/
Split & Split::operator=(
  const Split & other)  /**< is the Split object to be copied */
	{
    split_ntax		= other.split_ntax;
    nunits			= other.nunits;
    on_symbol		= other.on_symbol;
    off_symbol		= other.off_symbol;
    excl_symbol		= other.excl_symbol;
    mask			= other.mask;
	unit.resize(nunits, (split_t)0);
	std::copy(other.unit.begin(), other.unit.end(), unit.begin());
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing because all data members are automatically deleted.
*/
Split::~Split()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets non-constant data members using a pattern supplied in the form of the string `s' comprising a sequence of 
|   `on_symbol', `off_symbol' and (possibly) `excl_symbol' characters. The data members `on_symbol', `off_symbol' and
|   `excl_symbol' are not modified. For example, assuming `bits_per_unit' = 8, calling the function with the pattern 
|   "--*---***-" sets `split_ntax' to 10, `nunits' to 2,  `mask' to 3 (binary 00000011), `unit'[0] to 142 (binary 
|   10001110), and `unit'[1] to 0.
*/
void Split::CreateFromPattern(
  std::string s)    /**< is the string containing the pattern to use in constructing split */
    {
    unsigned slen = s.size();
    std::vector<unsigned> on_bits;
    for (unsigned k = 0; k < slen; ++k)
        {
        if (s[k] == on_symbol)
            on_bits.push_back(k);
        else if (s[k] == excl_symbol)
            excl_bits.push_back(k);
        else if (s[k] != off_symbol)
            throw XPhylogeny(str(boost::format("character in pattern (%c) not recognized as either the on symbol (%c), off symbol (%c) or excluded symbol (%c)") % s[k] % on_symbol % off_symbol % excl_symbol));
        }

    // No funny characters were found in the supplied pattern string, so we have a green light to build the split
    CalcNUnits(slen);   // sets split_ntax, nunits, and mask
	unit.resize(nunits, (split_t)0);
    for (std::vector<unsigned>::const_iterator it = on_bits.begin(); it != on_bits.end(); ++it)
        {
        unsigned k = *it;
        unsigned i = k/bits_per_unit;
        unsigned j = k % bits_per_unit;
        unit[i] |= (split_unity << j);
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes those data members (`split_ntax', `nunits', and `mask') that depend on the number of taxa in the tree. 
|   The data member `split_ntax' is set to the specified number of taxa `ntax'. The data member `nunits' equals the 
|   number of units (each holding `bits_per_unit' bits) needed to accommodate `split_ntax' taxa. Finally, low-order bits
|    in `mask' are set so that `mask' can be used to "see" only the relevant bits in the final unit (usually ntax does
|   not divide evenly into `bits_per_unit', so some bits are inevitably wasted and `mask' makes it easy to ensure that
|   we never pay attention to those bits.
*/
void Split::CalcNUnits(
  unsigned ntax)    /**< is the number of taxa in the tree */
	{
	if (ntax == 0)
		{
	    split_ntax = 0;
		nunits = 1;
		mask = (split_t)0;
		}
    else
        {
	    split_ntax = ntax;

        // Suppose bits_per_unit = 8
        //   7 taxa requires 1 unit  ==> 1 + (7-1)/8 = 1 + 6/8 = 1 + 0 = 1
        //   8 taxa requires 1 unit  ==> 1 + (8-1)/8 = 1 + 7/8 = 1 + 0 = 1
        //   9 taxa requires 2 units ==> 1 + (9-1)/8 = 1 + 8/8 = 1 + 1 = 2
	    nunits = 1 + ((split_ntax - 1)/bits_per_unit);

        // The remainder of this function is concerned with setting up the mask,
        // which has 0s for bits that are not used, and 1s for the ntax bits that
        // are used, in the final element of the units vector. For example, for 
        // bits_per_unit = 8, the last_bit value and the mask would look like this 
        // for several example ntax values:
        //
        //                     76543210
        //   ntax =  7  mask = 01111111  last_bit = 6
        //   ntax =  8  mask = 11111111  last_bit = 7
        //   ntax =  9  mask = 00000001  last_bit = 0
        //   ntax = 10  mask = 00000011  last_bit = 1
        //
	    unsigned num_unused_bits = (nunits*bits_per_unit) % split_ntax;
	    unsigned last_bit = bits_per_unit - num_unused_bits - 1;

        //@POL  I don't think this is necessary after all
	    // The following is necessary for the case in which split_ntax < bits_per_unit
	    //if (last_bit > ntax - 1)
		//    last_bit = ntax - 1;

	    mask = split_unity;
	    for (unsigned i = 0; i < last_bit; ++i) 
		    {
		    mask <<= 1;
		    mask |= split_unity;
		    }
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Serializes a Split object to a string `out', then returns a reference to `out'. Split objects output using this 
|   operator can be input using the corresponding operator>>. The information from `s' is space-delimited, and in the
|   following order: `bits_per_unit', `split_ntax', `on_symbol', `off_symbol', `excl_symbol', `unit[0]', `unit[1]', ..., 
|   `unit[nunits-1]'. Neither `nunits' nor `mask' are included because these can be recalculated from knowledge of
|   `bits_per_unit' and `split_ntax' using the function CalcNUnits. Note that `out' is cleared before any values are
|   appended.
*/
std::string & operator<<(
  std::string & out, /**< is the string on which to append the bit representation */
  const Split & s)   /**< is the Split object for which we want a string representation */
	{
    out = str(boost::format("%d %d %c %c %c") % s.bits_per_unit % s.split_ntax % s.on_symbol % s.off_symbol % s.excl_symbol);
	for (unsigned j = 0; j < s.nunits; ++j)
		{
        out += str(boost::format(" %d") % s.unit[j]);
		}
	return out;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inputs from supplied ifstream `in' all the information needed to fully construct the Split `s'. Use operator<< to 
|   output Split objects so that they can be read in using this operator. On failure, an XPhylogeny exception will be
|   thrown and `s' will remain unaffected.
*/
std::istream & operator>>(
  std::istream & in,    /**< is the ifstream from which to read the information needed for constructing a Split */
  Split & s)            /**< is the Split object that will be built from the information read from `in' */
	{
    Split tmp;
    unsigned b;
    in >> b >> tmp.split_ntax >> tmp.on_symbol >> tmp.off_symbol >> tmp.excl_symbol;

    // If any of the above read operations failed, bail out now
    if (!in)
        throw XPhylogeny("problem reading split dimensions from input stream");

    // Check to make sure the size of the type split_t is the same now as when the data was saved
    if (b != tmp.bits_per_unit)
        throw XPhylogeny(str(boost::format("the number of bits per unit differs on the system being used to read a split object (%d) compared to the system used to save the split object (%d)") % tmp.bits_per_unit % b));

    tmp.CalcNUnits(tmp.split_ntax); // builds mask and sets value of nunits
    tmp.unit.resize(tmp.nunits, (split_t)0);

	for (unsigned j = 0; j < s.nunits; ++j)
		{
		in >> tmp.unit[j];
		}

    if (!in)
        throw XPhylogeny("problem reading split data from input stream");

	return in;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a coded representation in which `on_symbol' represents bits that are "on" (i.e. set to 1), `off_symbol' 
|   represents "off" bits (i.e. set to 0), and `excl_symbol' represent bits that have been excluded. Uses function 
|   Split::CreateAndAppendPatternRepresentation to do the actual work.
*/
std::string Split::CreatePatternRepresentation() const
	{
	std::string s;
	s.reserve(split_ntax);
	CreateAndAppendPatternRepresentation(s);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a coded representation in which `on_symbol' represents bits that are "on" (i.e. set to 1), `off_symbol' 
|   represents "off" bits (i.e. set to 0), and `excl_symbol' represent bits that have been excluded. This function is
|   employed by Split::CreatePatternRepresentation. Note that taxon 0 occupies the leftmost position in the resulting
|   string, whereas in underlying representation things are not so tidy. Here is an example involving 10 taxa where 
|   `bits_per_unit' = 8 and taxa 1, 7 and 9 are "on":
|>
|   +--------+--------+
|   |10000010|00000010|
|   +--------+--------+
|     unit[0]  unit[1]
|>
|   This split would be output as follows, assuming `on_symbol' is an asterisk and `off_symbol' is a hyphen:
|>
|           -*-----*-*
|   taxon 0 ^        ^ taxon 9
|>
|   Note that "taxon 9" actually represents the 10th. taxon.
*/
void Split::CreateAndAppendPatternRepresentation(
  std::string & s) const	/**< is a reference to the string to which the representation will be appended */
	{
	if (split_ntax > 0)
		{
		unsigned ntax_added = 0;
		for (unsigned i = 0; i < nunits; ++i) 
			{
			for (unsigned j = 0; j < bits_per_unit; ++j) 
				{
                split_t bitmask = (split_unity << j);
				bool bit_is_set = ((unit[i] & bitmask) > (split_t)0);
                if (!excl_bits.empty() && std::binary_search(excl_bits.begin(), excl_bits.end(), ntax_added))
                    s += excl_symbol;
                else if (bit_is_set)
                    s += on_symbol;
                else
                    s += off_symbol;
				if (++ntax_added == split_ntax)
					return;
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a newick tree description representing the split and returns as a string. Each taxon is represented by an
|   integer starting at 1 (not 0!). The taxa representing "on" bits are listed first in the tree description. Thus, if
|   bits 1, 7 and 9 were "on", and 0, 2-6, and 8 were "off", the tree description would look like this:
|>
|   (2,8,10,(1,3,4,5,6,7,9))
|>
|   Note the lack of a terminating semicolon.
*/
std::string Split::CreateNewickRepresentation() const
	{
	std::string s;
	s += '(';
    // Spit out "on" bits first, but build a vector of off bits for later use
    std::vector<unsigned> off;
    unsigned ntax_added = 0;
    for (unsigned i = 0; i < nunits; ++i)
        {
		for (unsigned j = 0; j < bits_per_unit; ++j) 
			{
            split_t bitmask = (split_unity << j);
			bool bit_is_set = ((unit[i] & bitmask) > (split_t)0);
            //
            // i   j
            // 0   0  1  2  3  4  5  6  7
            // ==> 1  2  3  4  5  6  7  8 <== taxon_number
            //
            // i   j
            // 1   0  1  2  3  4  5  6  7
            // ==> 9 10 11 12 13 14 15 16 <== taxon_number
            //
            unsigned taxon_number = i*bits_per_unit + j + 1;
            if (bit_is_set)
                {
                s += str(boost::format("%d,") % taxon_number);
                }
            else
                off.push_back(taxon_number);
			if (++ntax_added == split_ntax)
				break;
			}
        }
	s += '(';

    // Now for the "off" bits
    std::vector<unsigned>::const_iterator it = off.begin();
    s += str(boost::format("%d") % (*it));
    for (++it; it != off.end(); ++it)
        {
        s += str(boost::format(",%d") % (*it));
        }
	s += "))";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a representation of the Split in which each unit is simply output as a decimal value. This representation is
|	not easily interpreted, as each unit comprises bits for several taxa (and it is the binary representation that makes 
|	sense). The units are output with unit[0] on the left. Here is an example involving bits_per_unit = 8:
|>
|   +--------+--------+
|   |10000010|00000010|
|   +--------+--------+
|     unit[0]  unit[1]
|>
|   The string produced from the above unit vector would be:
|>
|   130 2
|>
*/
std::string Split::CreateIdRepresentation() const
	{
	std::string s;
    PHYCAS_ASSERT(unit.size() > 0);
    s += str(boost::format("%d") % unit[0]);
	for (unsigned i = 1; i < nunits; ++i)
        s += str(boost::format(" %d") % unit[i]);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of bits that are set to one (i.e. are "on") in the Split.
*/
unsigned Split::CountOnBits() const
	{
    // This assumes that the unused bits in the last unit are all 0.
    unsigned num_bits_set = 0;

    for (SplitTVect::const_iterator i = unit.begin(); i != unit.end(); ++i)
        {
        // The  Kernighan method of counting bits is used below: see exercise 2-9 in the Kernighan and Ritchie book.
        // As an example, consider v = 10100010
        // c = 0:
        //   v     = 10100010
        //   v - 1 = 10100001
        //   ----------------
        //   new v = 10100000
        // c = 1:
        //   v     = 10100000
        //   v - 1 = 10011111
        //   ----------------
        //   new v = 10000000
        // c = 2:
        //   v     = 10000000
        //   v - 1 = 01111111
        //   ----------------
        //   new v = 00000000
        // c = 3:
        //   break out of loop because v = 0
        // 
        split_t v = *i;
        unsigned c = 0;
        for (; v; ++c)
            {
            v &= v - 1;
            }
        num_bits_set += c;
        }
    return num_bits_set;
    }
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of bits that are set to zero (i.e. are "off") in the Split. This function is slower than the
|   function CountOnBits unless `excl_bits' is empty, in which case the two functions will be essentially equal in speed.
|   The problem is that CountOnBits does not need to check for excluded bits because excluded bits are guaranteed to be
|   off. If the `excl_bits' vector is empty, then CountOffBits is slower only because it needs to invert a temporary
|   copy of each unit.
*/
unsigned Split::CountOffBits() const
	{
    // This assumes that the unused bits in the last unit are all 0.
    unsigned num_bits_unset = 0;

    // Pretend that there are no excluded bits
    for (SplitTVect::const_iterator i = unit.begin(); i != unit.end(); ++i)
        {
        // Set comment in CountOnBits about this Kernighan bit-counting algorithm
        split_t v = ~(*i);
        unsigned c = 0;
        for (; v; ++c)
            {
            v &= v - 1;
            }
        num_bits_unset += c;
        }

    if (!excl_bits.empty())
        {
        // Correct count to account for excluded bits
        for (std::vector<unsigned>::const_iterator i = excl_bits.begin(); i != excl_bits.end(); ++i)
            {
            // Subtract one for every bit that is unset (remember that bits that are unset in this split
            // are set in its inverse)
            if (!IsBitSet(*i))
                --num_bits_unset;
            }
        }

    // Account for irrelevant bits at the end
    unsigned num_irrelevant = (nunits*bits_per_unit - split_ntax);
    num_bits_unset -= num_irrelevant;

    return num_bits_unset;
    }
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns min(n,m), where n is the number of taxa on one side of the split and m is the number on the other side.
|	Trivial splits have m=1 or n=1, and thus have compexity 1, whereas the most complex split has complexity 
|	split_ntax/2 (note that this maximum holds whether or not split_ntax is even or odd).
*/
unsigned Split::CalcComplexity() const
	{
	const unsigned on_bits = CountOnBits();
	const unsigned off_bits = CountOffBits();
	return (on_bits < off_bits ? on_bits : off_bits);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Records the current values of the data members `bits_per_unit', `split_ntax', `nunits', and `mask', along with the
|   size of a split_t object on the current system, to a string object, which is returned.
*/
std::string Split::GetDimensionInfo()
	{
    return str(boost::format("Split information:\n  bits_per_unit = %d\n  split_ntax = %d\n  nunits = %d\n  mask = %d\n  sizeof(split_t) = %d") 
        % bits_per_unit % split_ntax % nunits % mask % sizeof(split_t));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if this->`unit' < `other'.`unit'.
*/
bool Split::IsLessThan(
  const Split & other) const /**< is the other Split object to which this Split object is being compared */
	{
    PHYCAS_ASSERT(unit.size() == other.unit.size());
    return (unit < other.unit);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if bit corresponding to taxon `t' is set, and returns false otherwise. The supplied value `t' should be
|   0 if it represents the first taxon.
*/
bool Split::IsBitSet(
	unsigned t) const /**< is the bit to consider, where `t' = 0 corresponds to first taxon */
	{
	PHYCAS_ASSERT(t < split_ntax);
	
	unsigned i = t/bits_per_unit;
	unsigned j = t % bits_per_unit;
	split_t x = unit[i] & (split_unity << j);
	return (x > (split_t)0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes `unit' vector and creates a new `unit' data member with `nunits' elements all initialized to 0. Assumes 
|   `nunits' is already set correctly, which will be the case if the member function Split::CalcNUnits has just been
|   called.
*/
void Split::Resize()
	{
	unit.resize(nunits, (split_t)0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the bit corresponding to taxon `t', where `t' = 0 corresponds to the first taxon.
*/
void Split::SetBit(
	unsigned t)	/**< the bit to set, where `t' = 0 corresponds to the first taxon */
	{
	PHYCAS_ASSERT(t < split_ntax);
	
	unsigned i = t/bits_per_unit;
	unsigned j = t % bits_per_unit;
	unit[i] |= (split_unity << j);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets all bits corresponding to values in the supplied vector `bits_to_set'. Values in `bits_to_set' are 0-based
|   indices (i.e. use 0 to set the first bit).
*/
void Split::SetBits(
  const std::vector<unsigned> bits_to_set)	/**< the bits to set */
	{
    for (std::vector<unsigned>::const_iterator it = bits_to_set.begin(); it != bits_to_set.end(); ++it)
        {
        unsigned value = *it;
	    PHYCAS_ASSERT(value < split_ntax);
    	unsigned i = value/bits_per_unit;
	    unsigned j = value % bits_per_unit;
	    unit[i] |= (split_unity << j);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Unsets the bit corresponding to taxon `t', where `t' = 0 corresponds to the first taxon.
*/
void Split::UnsetBit(
  unsigned t)	/**< the bit to set, where `t' = 0 corresponds to the first taxon */
	{
	PHYCAS_ASSERT(t < split_ntax);
	unsigned i = t/bits_per_unit;	 
	unsigned j = t % bits_per_unit; 
	unit[i] &= (~(split_unity << j));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Unsets (i.e. clears) all bits corresponding to values in the supplied vector `bits_to_unset'. Values in 
|   `bits_to_unset' are 0-based indices (i.e. use 0 to unset the first bit).
*/
void Split::UnsetBits(
  const std::vector<unsigned> bits_to_unset)	/**< the bits to unset */
	{
    for (std::vector<unsigned>::const_iterator it = bits_to_unset.begin(); it != bits_to_unset.end(); ++it)
        {
        unsigned value = *it;
	    PHYCAS_ASSERT(value < split_ntax);
    	unsigned i = value/bits_per_unit;
	    unsigned j = value % bits_per_unit;
    	unit[i] &= (~(split_unity << j));
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If this split is subsumed in `other', then a bitwise AND for any unit should equal the unit from this split. If this
|	is not the case, it means there was a 0 in other for a bit that is set in this split. For example:
|>
|	**-*--**-*-  <-- this split
|	*****-****-  <-- other split
|	----------------------------
|	**-*--**-*-  <-- bitwise AND
|>
|	In the above example, this split is subsumed in the other split because the bitwise AND equals this split. 
|	Note that every split is by definition subsumed in itself. The `start_unit' argument is used by Split::IsCompatible.
*/
bool Split::SubsumedIn(
	const Split & other,			/**< the split for comparison */
	unsigned start_unit) const	    /**< elements before unit[start_unit] will not be considered */
	{
	for (unsigned i = start_unit; i < nunits; ++i)
		{
		if ((unit[i] & other.unit[i]) != unit[i]) 
			return false;
		}
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if this split and `other' are compatible. The two splits a and b are compatible if a & b is nonzero
|   and also not equal to either a or b. For example, these two splits 
|>
|   split a: -***---*--
|   split b: ----***--*
|     a & b: ----------
|>
|   are compatible, because a & b = 0. The two splits below are also compatible because a & b == b:
|>
|   split a: -****-*---
|   split b: --**--*--- <-
|     a & b: --**--*--- <-
|>
|   These two splits, on the other hand, are not compatible because a & b != 0 and is not equal to either a or b:
|>
|   split a: -***---*--
|   split b: ---***---*
|     a & b: ---*------ 
|>
*/
bool Split::IsCompatible(
  const Split & other) const	/**< the split for comparison */
	{
	for (unsigned i = 0; i < nunits; ++i)
		{
        split_t a       = unit[i];
        split_t b       = other.unit[i];
		split_t a_and_b = (a & b);
        bool equals_a   = (a_and_b == a);
        bool equals_b   = (a_and_b == b);
        if (a_and_b && !(equals_a || equals_b))
            {
            // A failure of any unit to be compatible makes the entire split incompatible
            return false;
            }

        //@POL not sure what the following was all about; perhaps if the new code above fails to work
        // I'll find that the original code below was actually correct! At the moment, it makes no sense
        // to me at all, however.
        //
		//if (t != unit[i])
		//	{
		//	if (t != other.unit[i])
        //        {
		//		return false;
        //        }
		//	return (++i == nunits ? true : other.SubsumedIn(*this, i));
		//	}
		//else if (t != other.unit[i]) 
		//	return (++i == nunits ? true : SubsumedIn(other, i));	
		}

    // None of the units were incompatible, so that means the splits are compatible
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts each element of the unit array to its bitwise complement, then clears the irrelevant bits in the final 
|   element. The irrelevant bits arise because bits must be allocated in chunks of size `bits_per_unit'. For example, 
|   suppose `split_ntax' = 5, `bits_per_unit' = 4, and thus `nunits' = 2. The `unit' vector is shown below for a split
|   in which the bits for the last three taxa are "on":
|>
|	+---+---+---+---+  +---+---+---+---+
|	| 0 | 0 | 0 | 1 |  | 1 | 1 | 0 | 0 |
|	+---+---+---+---+  +---+---+---+---+
|	      unit[1]           unit[0]
|>
|	The three most significant (i.e. left-most) bits in unit[1] in this case are irrelevant because there are 8 bits
|	total and only 5 of them are used to store information. The mask is used to clear these irrelevant bits in unit[1].
|	In this example, mask would look like this:
|>
|	+---+---+---+---+
|	| 0 | 0 | 0 | 1 |
|	+---+---+---+---+
|	      mask
|>
|	Note that ordinarily `bits_per_unit' would be greater than 4 (32 or 64 are common values for 32-bit or 64-bit
|   systems, respectively). The `unit' vector after this method is called would be:
|>
|	+---+---+---+---+  +---+---+---+---+
|	| 0 | 0 | 0 | 0 |  | 0 | 0 | 1 | 1 |
|	+---+---+---+---+  +---+---+---+---+
|	      unit[1]           unit[0]
|>
*/
void Split::InvertSplit()
	{
    for (SplitTVect::iterator i = unit.begin(); i != unit.end(); ++i) 
		{
		split_t x = *i;
		*i = ~x;
		}
	unit[nunits - 1] &= mask;

    // Must ensure that none of the bits now set corresponds to an excluded bit (cannot use the fast version of
    // CountOnBits unless all excluded bits have been cleared)
    if (!excl_bits.empty())
        UnsetBits(excl_bits);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	For each split in the `ref_tree', checks to see if the split is in the `test_tree' SplitSet. If not, inserts the
|	split into `missing'. Assumes `ref_tree' and `test_tree' SplitSet objects are already filled.
*/
unsigned FindSplitsAbsentInTestTree(
	const SplitSet & ref_tree,	/**< is the set of splits to check */
	const SplitSet & test_tree,	/**< is the set of splits to test */
	SplitSet & missing)			/**< in the end, contains splits present in `ref_tree' but absent from `test_tree' */
	{
    missing.clear();
	SplitSet::iterator insertLoc = missing.begin();
	unsigned num_absent_splits = 0;
	for (SplitSet::const_iterator it = ref_tree.begin(); it != ref_tree.end(); ++it)
		{
		if (test_tree.find(*it) == test_tree.end())
			{
			insertLoc = missing.insert(insertLoc, *it);
			++num_absent_splits;
			}
		}
	return num_absent_splits;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns true if and only if this split is not equal to `other'.
*/
bool Split::operator!=(
  const Split & other) const    /**< is the split for comparison */
    {
    if (*this == other)
        return false;
    else
        return true;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a split in which each element of the `unit' vector is its bitwise AND with the corresponding element from
|   the `unit' vector of `other'. For example:
|>
|   split a: -*-**-***
|   split b: -****----
|     a & b: -*-**----
|>
*/
Split Split::operator&(
  const Split & other) const    /**< is the split for comparison */
	{
	Split tmp(*this);
	return tmp &= other;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|   Returns a reference to this Split object after replacing each element of its `unit' vector with the bitwise AND of
|   that element with the corresponding element of the unit vector of `other'. For example:
|>
|   split a: -*-**-***
|   split b: -****----
|         a: -*-**---- after a &= b
|>
*/
Split & Split::operator&=(
  const Split & other)   /**< is the split for comparison */
	{
	for (unsigned i = 0; i < nunits; ++i)
		unit[i] &= other.unit[i];
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a split in which each element of the `unit' vector is its bitwise XOR with the corresponding element from
|   the `unit' vector of `other'. For example:
|>
|   split a: -*-**-***
|   split b: -****----
|     a ^ b: --*---***
|>
*/
Split Split::operator^(
  const Split & other) const   /**< is the split for comparison */
    {
    Split tmp(*this);
    return tmp ^= other;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a reference to this Split object after replacing each element of its `unit' vector with the bitwise XOR of
|   that element with the corresponding element of the unit vector of `other'. For example:
|>
|   split a: -*-**-***
|   split b: -****----
|         a: --*---*** after a ^= b
|>
*/
Split & Split::operator^=(
  const Split & other)   /**< is the split for comparison */
    {
    for (unsigned i = 0; i < nunits; ++i)
        unit[i] ^= other.unit[i];
    return *this;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a split in which each element of the `unit' vector is its bitwise OR with the corresponding element from
|   the `unit' vector of `other'. For example:
|>
|   split a: -*-**-***
|   split b: -****----
|     a | b: -****-***
|>
*/
Split Split::operator|(
  const Split & other) const   /**< is the split for comparison */
	{
	Split tmp(*this);
	return tmp |= other;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a reference to this Split object after replacing each element of its `unit' vector with the bitwise OR of
|   that element with the corresponding element of the unit vector of `other'. For example:
|>
|   split a: -*-**-***
|   split b: -****----
|         a: -****-*** after a |= b
|>
*/
Split & Split::operator|=(
  const Split & other)   /**< is the split for comparison */
	{
	for (unsigned i = 0; i < nunits; ++i)
		unit[i] |= other.unit[i];
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls IsLessThan(`other') and returns the value returned by that function.
*/
bool Split::operator<(
  const Split & other) const   /**< is the split for comparison */
	{
	return IsLessThan(other);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Equals(`other') and returns the value returned by that function.
*/
bool Split::operator==(
  const Split & other) const   /**< is the split for comparison */
	{
	return Equals(other);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if all elements of the unit vector are equal to all elements of the Split object `other'.
*/
bool Split::Equals(
  const Split & other) const   /**< is the split for comparison */
	{
	return std::equal(unit.begin(), unit.end(), other.unit.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a bitwise OR of each element of unit with the corresponding element in the unit vector of `other'. The
|   result is a union of the sets defined by the two Split objects, and is useful in creating the Split for an interior
|   node, which is the union of the splits of its immediate descendants.
*/
void Split::CombineWith(
  const Split & other)   /**< is the split to combine with this Split object */
	{
	*this |= other;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts this split into the intersection between this split and the other split. The other split remains unaffected
|	but this Split object is modified such that each element in the unit array undergoes a bitwise AND operation with
|	the corresponding element from unit vector of `other'.
*/
void Split::IntersectWith(
	const Split & other)	/* the split with which to intersect */
	{
	*this &= other;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `split_ntax'.
*/
unsigned Split::GetNTaxa() 
    {
    return split_ntax;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns vector of unsigned integers each of which represents a bit that is set in this split (first is 0). This
|   function removes values corresponding to excluded taxa before returning the list.
*/
std::vector<unsigned> Split::GetOnList() const
    {
    unsigned k = 0;
    std::vector<unsigned> v;
    for (unsigned i = 0; i < nunits; ++i)
        {
        for (unsigned j = 0; j < bits_per_unit; ++j)
            {
            split_t bit = (split_unity << j);
            bool is_on = ((unit[i] & bit) > (split_t)0);
            if (is_on)
                v.push_back(k);
            if (++k == split_ntax)
                break;
            }
        }

    // Eliminate values also in excl_bits vector
    if (!excl_bits.empty())
        {
        std::vector<unsigned>::iterator last = std::set_difference(v.begin(), v.end(), excl_bits.begin(), excl_bits.end(), v.begin());
        v.erase(last, v.end());
        }
    return v;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns vector of unsigned integers each of which represents a bit that is not set in this split (first is 0). This
|   function removes values corresponding to excluded taxa before returning the list.
*/
std::vector<unsigned> Split::GetOffList() const
    {
    unsigned k = 0;
    std::vector<unsigned> v;
    for (unsigned i = 0; i < nunits; ++i)
        {
        for (unsigned j = 0; j < bits_per_unit; ++j)
            {
            split_t bit = (split_unity << j);
            bool is_on = ((unit[i] & bit) > (split_t)0);
            if (!is_on)
                v.push_back(k);
            if (++k == split_ntax)
                break;
            }
        }

    // Eliminate values also in excl_bits vector
    if (!excl_bits.empty())
        {
        std::vector<unsigned>::iterator last = std::set_difference(v.begin(), v.end(), excl_bits.begin(), excl_bits.end(), v.begin());
        v.erase(last, v.end());
        }
    return v;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns vector of unsigned integers each of which represents a bit that is currently excluded (first is 0). The bits
|   that are excluded are not considered either on or off, and are stored in the data member `excl_bits'.
*/
std::vector<unsigned> Split::GetExcludedList() const
    {
    unsigned k = 0;
    std::vector<unsigned> v;
    if (excl_bits.empty())
        return v;
    else
        {
        v.resize(excl_bits.size(), 0);
        std::copy(excl_bits.begin(), excl_bits.end(), v.begin());
        }
    return v;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Establishes a list of bits that are currently excluded. These bits will always be off, but they should not be 
|   returned in the list of values produced by GetOffList. Note that this function deletes any existing elements in the
|   `excl_bits' vector.
*/
void Split::SetExcluded(const std::vector<unsigned> excl)
    {
    if (excl.empty())
        {
        excl_bits.clear();
        }
    else
        {
        excl_bits.resize(excl.size(), 0);
        std::copy(excl.begin(), excl.end(), excl_bits.begin());

        // Must ensure excl_bits is sorted because binary_search algorithm (used in CreateAndAppendPatternRepresentation)
        // and set_difference algorithm (used in GetOnBits and GetOffBits) require it
        std::sort(excl_bits.begin(), excl_bits.end()); 

        // Also ensure that none of the bits now set corresponds to an excluded bit (cannot use the fast version of
        // CountOnBits unless all excluded bits have been cleared)
        UnsetBits(excl_bits);
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Set the value of the data member `on_symbol' used to represent set bits by functions such as 
|   Split::CreatePatternRepresentation.
*/
void Split::SetOnSymbol(const char c)
    {
    on_symbol = c;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Set the value of the data member `off_symbol' used to represent unset bits by functions such as 
|   Split::CreatePatternRepresentation.
*/
void Split::SetOffSymbol(const char c)
    {
    off_symbol = c;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Set the value of the data member `excl_symbol' used to represent excluded bits by functions such as 
|   Split::CreatePatternRepresentation.
*/
void Split::SetExcludedSymbol(const char c)
    {
    excl_symbol = c;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns current value of the data member `on_symbol' used to represent set bits by functions such as 
|   Split::CreatePatternRepresentation.
*/
char Split::GetOnSymbol() const
    {
    return on_symbol;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns current value of the data member `off_symbol' used to represent unset bits by functions such as 
|   Split::CreatePatternRepresentation.
*/
char Split::GetOffSymbol() const
    {
    return off_symbol;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns current value of the data member `excl_symbol' used to represent excluded bits by functions such as 
|   Split::CreatePatternRepresentation.
*/
char Split::GetExcludedSymbol() const
    {
    return excl_symbol;
    }
