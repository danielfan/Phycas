#include "phycas/force_include.h"
#include <cmath>
#include "ncl/misc/nxs_index_set.hpp"
#include "ncl/taxa/nxs_taxa_manager.hpp"
#include "phycas/trees/split.hpp"
#if defined(POL_PHYCAS)
#	include "pyphy/prob_dist/basic_cdf.hpp"
#	include "pyphy/prob_dist/basic_lot.hpp"
#else
#	include "phycas/rand/lot.hpp"
#endif
using std::string;

//const unsigned Split::bits_per_unit = (CHAR_BIT)*sizeof(split_t);	/* number of bits in a variable of type split_t */
//const split_t  Split::unity			= (split_t)1;
	
unsigned 	Split::splitNTax		= 4; // the number of taxa that the split represents

unsigned	Split::nunits			= 1;
char		Split::on_symbol		= '*';
char		Split::off_symbol		= '-';

	split_t		Split::mask				= (split_t)0;

#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::log;
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes those STATIC_DATA data members (splitNTax, nunits) that depend on the number of taxa. splitNTax is set to the 
|	specified number of taxa, n. nunits equals the number of units (each holding bits_per_unit bits) needed to 
|	accommodate splitNTax taxa.
*/
void Split::CalcStatics(unsigned n)
	{
	if (n == 0)
		{
	    Split::splitNTax = 0;
		Split::nunits = 1;
		Split::mask = (split_t)0;
		return;
		}

	Split::splitNTax = n;
	Split::nunits = 1 + ((Split::splitNTax - 1) / Split::bits_per_unit);

	unsigned num_irrelevant_bits = (Split::nunits*Split::bits_per_unit) % Split::splitNTax;
	unsigned last_bit = Split::bits_per_unit - num_irrelevant_bits - 1;

	// The following is necessary for the case in which splitNTax < bits_per_unit
	//
	if (last_bit > n - 1)
		last_bit = n - 1;

	Split::mask = (split_t)1;
	for (unsigned i = 0; i < last_bit; ++i) 
		{
		mask <<= 1;
		mask |= Split::unity;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Outputs all nunits values of the unit array of s to the string out, with each value preceded by a blank space, 
|	then returns a reference to out. Split objects output using this operator can be input using the corresponding >> 
|	operator.
*/
string &operator<<(string &out, const Split &s)
	{
	//@POL not sure how useful this is because it doesn't output all the information about the split that might be needed
	for (unsigned j = 0; j < Split::nunits; ++j)
		{
		out << ' ' << s.unit[j];	
		}
	return out;
	}

void Split::SetBits(const NxsIndexSet & taxset)
	{
	for (NxsIndexSet::const_iterator iter = taxset.begin(); iter != taxset.end(); ++iter)
			SetBit(*iter);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates unit array nunits long and calls Clear().
*/
Split::Split(const NxsIndexSet &s)
	{
	unit = new split_t[nunits];
	Clear();
	SetBits(s);
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Inputs from supplied ifstream all nunits values of the unit array of s. Use operator<<to output Split objects 
|	so that they can be	read in using this operator.
*/
std::istream &operator>>(std::istream &in, Split &s)
	{
	//@POL not sure how useful this is because it doesn't input all the information about the split that might be needed
	for (unsigned j = 0; j < Split::nunits; ++j)
		{
		in >> s.unit[j];
		}
	return in;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a coded representation in which Split::on_symbol represents on bits and Split::off_symbol represents off 
|	bits. Uses Split::CreateAndAppendPatternRepresentation() to do the actual work.
*/
string Split::CreatePatternRepresentation() const
	{
	string s;
	s.reserve(Split::splitNTax);
	CreateAndAppendPatternRepresentation(&s);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a newick tree description representing the split and returns as a string.
*/
string Split::CreateNewickRepresentation(const NxsTaxaManager *taxaMgr) const
	{
	NxsIndexSet includedTaxa, excludedTaxa;
	GetIncludedExcludedRepresentation(&includedTaxa, &excludedTaxa);

	string s;
	if (includedTaxa.size() < 2 || excludedTaxa.size() < 2) 
		{
		s << '(';
		bool first = true;
		if (excludedTaxa.size() > 0)
			{
			for (NxsIndexSet::const_iterator exclIt = excludedTaxa.begin(); exclIt != excludedTaxa.end(); ++exclIt)
				{
				string n = taxaMgr->GetLabel(*exclIt);
				if (!first)
					s << ',';
				s << ConvertToNexusToken(n); 
				first = false;
				}
			}
		if (includedTaxa.size() > 0)
			{
			for (NxsIndexSet::const_iterator inclIt = includedTaxa.begin(); inclIt != includedTaxa.end(); ++inclIt)
				{
				string n = taxaMgr->GetLabel(*inclIt);
				if (!first)
					s << ',';
				s << ConvertToNexusToken(n); 
				first = false;
				}
			}
		s << ')';
		}
	else	// not a trivial split
		{
		s << '(';
		for (NxsIndexSet::const_iterator exclIt = excludedTaxa.begin(); exclIt != excludedTaxa.end(); ++exclIt)
			{
			string n = taxaMgr->GetLabel(*exclIt);
			s << ConvertToNexusToken(n) << ',';
			}
		s << '(';
		bool first = true;
		for (NxsIndexSet::const_iterator inclIt = includedTaxa.begin(); inclIt != includedTaxa.end(); ++inclIt)
			{
			string n = taxaMgr->GetLabel(*inclIt);
			if (!first)
				s << ',';
			s << ConvertToNexusToken(n); 
			first = false;
			}
		s << "))";
		}

	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts 0-offset index of taxon into NxsIndexSet inc if bit for this taxon is currently set; otherwise, inserts
|	taxon index into NxsIndexSet exc. Note: clears inc and exc before starting.
*/
void Split::GetIncludedExcludedRepresentation(NxsIndexSet *inc, NxsIndexSet *exc) const
	{
	inc->clear();
	exc->clear();
	for (unsigned i = 0; i < splitNTax; ++i) 
		{
		if (IsBitSet(i))
			inc->insert(i);
		else
			exc->insert(i);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a representation of the split in which each unit is simply output in decimial. This representation is not
|	easily interpreted, as each unit comprises bits for several taxa (and it is the binary representation that makes 
|	sense). Note: the units are output with unit[0] on the left.
*/
string Split::CreateIdRepresentation() const
	{
	//@POL Mark, I changed i = 1 to i = 0 (I don't know why my original code started at 1)
	// and added a space between each number

	string s;
	for (unsigned i = 0; i < nunits; ++i)
		s << ' ' << unit[i];
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of bits that are set to one in the split representation
*/
unsigned Split::CountOnBits() const
	{
	//@could be sped up using on bit lookup
	//@POL this counts irrelevant bits (i.e. those not representing any taxon) - shouldn't we mask off these nonsensical bits
	unsigned nOnBits = 0;
	for (unsigned i = 0; i < nunits; ++i) 
		{
		split_t k = Split::unity;
		for (unsigned j = 0; j < Split::bits_per_unit; ++j)
			{
			if ((unit[i] & k) != (split_t) 0)
				++nOnBits;
			k <<= 1;
			}
		}
	return nOnBits;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns min(n,m), where n is the number of taxa on one side of the split and m is the number on the other side.
|	Trivial splits have m=1 or n=1, and thus have compexity 1, whereas the most complex split has complexity 
|	splitNTax/2 (note that this maximum holds whether or not splitNTax is even or odd).
*/
unsigned Split::CalcComplexity() const
	{
	const unsigned complexity = CountOnBits();
	const unsigned complementComplexity = Split::splitNTax - complexity;
	return (complexity > complementComplexity ? complementComplexity : complexity);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Creates a coded representation in which the character Split::on_symbol represents on bits and Split::off_symbol
|	represents off bits. This function is employed by Split::CreatePatternRepresentation. Note that the order in which
|	bits are presented is taxon 0 leftmost, whereas in underlying representation things are not so tidy (units are laid
|	out in memory left to right, but within unit[0], for instance, least significant bit represents taxon 0).
*/
void Split::CreateAndAppendPatternRepresentation(
	string *s) const	/* pointer to string object to which the representation will be appended */
	{
	assert(s != NULL); //@pol -> mth: why not use string & here?
	if (Split::splitNTax > 0)
		{
		unsigned nTaxAdded = 0;
		for (unsigned i = 0;; ++i) 
			{
			for (unsigned j = 0; j < Split::bits_per_unit; ++j) 
				{
				bool bit_is_set = ((unit[i] & (Split::unity << j)) > (split_t)0);
				*s << (bit_is_set ? Split::on_symbol : Split::off_symbol);
				if (++nTaxAdded == Split::splitNTax)
					return;
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Split prior probability given a flat prior on tree topologies equals:
|>
|     (no. rooted trees with n tips) * (no. rooted trees with m tips)
|     ---------------------------------------------------------------
|                  (no. unrooted trees with n + m tips)
|
|     [ (2n-3)! / 2^(n-2) (n-2)! ] [ (2m-3)! / 2^(m-2) (m-2)! ]
|     ---------------------------------------------------------
|              [ (2n + 2m - 5)! / 2^(n+m-3) (n+m-3)!]
|>
|	where n is the number of taxa on one side of split and m is the number of taxa on the other side of the split.
|	This function returns the natural logarithm of the split prior probability for this split.
*/
double Split::CalcLnPriorProb(unsigned n, unsigned m)
	{
	if (n < 2 || m < 2)
		return 0.0;

	double ln_prior = log(2.0);

#if defined(POL_PHYCAS)
	phycas::CDF cdf;
	ln_prior += cdf.LnGamma((double)(2*n - 2));
	ln_prior += cdf.LnGamma((double)(2*m - 2));
	ln_prior += cdf.LnGamma((double)(n + m - 2));
	ln_prior -= cdf.LnGamma((double)(2*(n + m) - 4));
	ln_prior -= cdf.LnGamma((double)(m - 1));
	ln_prior -= cdf.LnGamma((double)(n - 1));
#else
	ln_prior += Lot::LotLnGamma((double)(2*n - 2));
	ln_prior += Lot::LotLnGamma((double)(2*m - 2));
	ln_prior += Lot::LotLnGamma((double)(n + m - 2));
	ln_prior -= Lot::LotLnGamma((double)(2*(n + m) - 4));
	ln_prior -= Lot::LotLnGamma((double)(m - 1));
	ln_prior -= Lot::LotLnGamma((double)(n - 1));
#endif
	return ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Records the current values of the STATIC_DATA data members to a string object, which is returned. Useful for debugging.
*/
string Split::GetStaticsInfo()
	{
	string s = "Split statics information:";
	s << "\n  Split::bits_per_unit = ";
	s << Split::bits_per_unit;
	s << "\n  Split::splitNTax     = ";
	s << Split::splitNTax;
	s << "\n  Split::nunits        = ";
	s << Split::nunits;
	s << "\n  Split::mask          = ";
	s << (long)Split::mask;
	s << "\n  sizeof(split_t)      = ";
	s << (long)sizeof(split_t);
	s << "\n  sizeof(double)       = ";
	s << (long)sizeof(double);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if this split is less than other, and false otherwise. If there is only one element in the unit array,
|	this split is less than the other split if unit[0] < other.unit[0]. If the unit array has size > 1, then elements
|	are visited starting with unit[nunits - 1], and this split is declared less than the other split if any element
|	is discovered to be less than the corresponding element in other's unit array. False is returned if any element is
|	discovered to be greater than the corresponding element in other's unit array, or if all elements are equal.
*/
bool Split::IsLessThan( const Split& other ) const
	{
	//@pol what about irrelevant bits? 
	for (unsigned i = nunits - 1;; --i) 
		{
		if (unit[i] < other.unit[i]) 
			return true;
		if (unit[i] > other.unit[i] || i == 0)
			return false;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if bit corresponding to taxon t is set, and returns false otherwise.
*/
bool Split::IsBitSet(
	unsigned t) const /* the 0-offset bit to consider, where t = 0 corresponds to first taxon */
	{
	NXS_ASSERT(unit != NULL);
	NXS_ASSERT(t < Split::splitNTax);
	
	const unsigned i = t/bits_per_unit;
	const unsigned j = t % bits_per_unit;
	const split_t x = unit[i] & (Split::unity << j);
	return (x > (split_t)0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes unit array and creates a new unit array with nunits elements. Assumes nunits is already set correctly.
*/
void Split::Resize()
	{
	delete [] unit;
	unit = new split_t[nunits];
	}

//BUGNOTE: POL 1/1/2002
//	Split::unity STATIC_DATA variable created because (1 << 32) yielded 0 on 64-bit machine
//	owing to the fact that the 1 was a 32-bit int rather than a 64-bit long. Could have
//	fixed by saying (1L << 32) but thought it was cleaner to declare a variable of 
//	exactly the right type (i.e. split_t) to use in these situations.
//
// t is taxon number starting with 0 for first taxon

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the bit corresponding to taxon t.
*/
void Split::SetBit(
	unsigned t)	/* the 0-offset bit to set, where t = 0 corresponds to the first taxon */
	{
	NXS_ASSERT(unit != NULL);
	NXS_ASSERT(t < Split::splitNTax);
	
	unsigned i = t/bits_per_unit;	//@POL t >> (log2_bits_per_unit); 
	unsigned j = t % bits_per_unit; //@POL t & (bits_per_unit - 1)
	unit[i] |= (Split::unity << j);
	}



//DANIEL 2-Dec-2004
/*----------------------------------------------------------------------------------------------------------------------
|	Unsets the bit corresponding to taxon t.
*/
void Split::UnsetBit(
	unsigned t)	/* the 0-offset bit to set, where t = 0 corresponds to the first taxon */
	{
	NXS_ASSERT(unit != NULL);
	NXS_ASSERT(t < Split::splitNTax);
	unsigned i = t/bits_per_unit;	 
	unsigned j = t % bits_per_unit; 
	unit[i] &= (~(Split::unity << j));
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
|	Note that every split is by definition subsumed in itself. The startUnit argument is used by Split::IsCompatible.
*/
bool Split::SubsumedIn(
	const Split &other,			/* the split for comparison */
	unsigned startUnit) const	/* elements before unit[startUnit] will not be considered */
	{
	for (unsigned i = startUnit; i < nunits; ++i)
		{
		if ((unit[i] & other.unit[i]) != unit[i]) 
			return false;
		}
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if this split and the argument are compatible, which means this split is subsumed in other and other is
|	subsumed in this split. See Split::SubsumedIn for definition of "subsumed in".
*/
bool Split::IsCompatible(
	const Split &other) const	/* the split for comparison */
	{
	for (unsigned i = 0; i < nunits; ++i)
		{
		const split_t t = (unit[i] & other.unit[i]);
		if (t != unit[i])
			{
			if (t != other.unit[i]) 
				return false;
			return (++i == nunits ? true : other.SubsumedIn(*this, i));
			}
		else if (t != other.unit[i]) 
			return (++i == nunits ? true : SubsumedIn(other, i));	
		}
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts each element of the unit array to its bitwise complement, then clears the irrelevant bits in 
|	the final element. The irrelevant bits arise because bits must be allocated in chunks of size bits_per_unit. For
|	example, suppose splitNTax = 5, bits_[er_unit = 4, and thus nunits = 2. The unit array is shown below:
|>
|	+---+---+---+---+  +---+---+---+---+
|	| - | - | - | 4 |  | 3 | 2 | 1 | 0 |
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
|	Note that ordinarily bits_per_unit would be greater than 4 (32 or 64 are more likely values).
*/
void Split::InvertSplit()
	{
	// Note: must use int here (or else choose another termination condition)
	//
	for (int i = (int)nunits - 1; i >= 0; --i) 
		{
		split_t x = unit[i];
		unit[i] = ~x;
		}
	unit[nunits - 1] &= Split::mask;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	For each split in the refTree SplitSet, checks to see if the split is in the testTree SplitSet. If not, inserts the
|	split into the missing SplitSet. Assumes refTree and testTree SplitSet objects are already filled.
*/
unsigned FindSplitsAbsentInTestTree(
	const SplitSet &refTree,	/* see description for missing */
	const SplitSet &testTree,	/* see description for missing */
	SplitSet *missing)			/* in the end, contains splits present in refTree SplitSet but absent from testTree SplitSet */
	{
	//pol -> mth: check documentation for accuracy
	//pol -> mth: shouldn't missing be cleared before starting?
	SplitSet::iterator insertLoc = missing->begin();
	unsigned nAbsentSplits = 0;
	for (SplitSet::const_iterator rtIt = refTree.begin(); rtIt != refTree.end(); ++rtIt)
		{
		if (testTree.find(*rtIt) == testTree.end())
			{
			insertLoc = missing->insert(insertLoc, *rtIt);
			++nAbsentSplits;
			}
		}
	return nAbsentSplits;
	}

#if 0
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
bool Split::Compatible( const Split& other ) const
	{
	// true if this->SubsumedIn(other) || other.SubsumedIn(*this)
	for (unsigned i = 0; i < nunits; ++i)
		{
		if (unit[i] != other.unit[i])
			{
			// we've reached the first element that is not equal in this and the other see if one is subsumed in the other.
			//
			const split_t x = (unit[i] & other.unit[i]);
			if (x == unit[i])
				{
				for (i++; i < nunits; ++i)
					{
					const split_t y = (unit[i] & other.unit[i]);
					if (x != unit[i])
						return false;
					}
				return true;
				}
			if (x == other.unit[i])
				{
				for (i++; i < nunits; ++i)
					{
					const split_t z = (unit[i] & other.unit[i]);
					if (x != other.unit[i])
						return false;
					}
				return true;
				}
			return false;
			}
		}
	return true;	//they are equal
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void Split::CreatePatternRepresentation( char* s, unsigned slen ) const
	{
	assert( slen > Split::ntax );
	
	unsigned pos = 0;
	for( unsigned i = 0; i < nunits && pos < Split::ntax; ++i ) 
		{
		split_t x = unit[i];
		for( unsigned j = 0; j < Split::bits_per_unit; ++j ) 
			{
			bool is_set = ( ( x & (Split::unity << j) ) > 0L );
			if( is_set )
				s[pos++] = '*';
			else
				s[pos++] = '-';
			if( pos >= Split::ntax )
				break;
			}
		}
	
	s[pos] = '\0';
	}


#endif
