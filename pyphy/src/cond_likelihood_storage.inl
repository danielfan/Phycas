#if ! defined(COND_LIKELIHOOD_STORAGE_INL)
#define COND_LIKELIHOOD_STORAGE_INL

//#define OBSOLETE_DEBUGGING_CODE
#include "pyphy/src/xlikelihood.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `realloc_min' to 1 and `num_patterns', `num_rates' and `num_states' all to 0. 
*/
inline CondLikelihoodStorage::CondLikelihoodStorage()
  : 
  num_patterns(0), 
  num_rates(0), 
  num_states(0), 
  num_created(0),
  realloc_min(1)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor calls clearStack to destroy all CondLikelihoodShPtr objects currently stored.
*/
inline CondLikelihoodStorage::~CondLikelihoodStorage()
	{
	clearStack();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `num_patterns', which holds the number of patterns used in determining the size
|	of newly-created CondLikelihood objects.
*/
inline unsigned CondLikelihoodStorage::getNumPatterns() const
	{
	return num_patterns;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `num_rates', which holds the number of rates used in determining the size of
|	newly-created CondLikelihood objects.
*/
inline unsigned CondLikelihoodStorage::getNumRates() const
	{
	return num_rates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `num_states', which holds the number of states used in determining the size of
|	newly-created CondLikelihood objects.
*/
inline unsigned CondLikelihoodStorage::getNumStates() const
	{
	return num_states;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of bytes allocated for each CLA. This equals sizeof(LikeFltType) times the product of the number of
|	patterns, number of rates and number of states.
*/
inline unsigned CondLikelihoodStorage::bytesPerCLA() const
	{
	return (unsigned)(num_patterns*num_rates*num_states*sizeof(LikeFltType));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of the data member `num_created', which keeps track of the number of CondLikelihood objects created 
|	since this object was constructed, or since the last call to clearStack, which resets `num_created' to zero.
*/
inline unsigned CondLikelihoodStorage::numCLAsCreated() const
	{
	return num_created;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current number of elements stored in the data member `cl_stack'. The total number of CLAs currently
|	checked out to the tree can be obtained as `num_created' minus the value returned by this function.
*/
inline unsigned CondLikelihoodStorage::numCLAsStored() const
	{
	return (unsigned)cl_stack.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `realloc_min' to 1, `cl_len' to `cond_like_len', and then calls the CondLikelihoodStorage::fillTo 
|	member function to add `starting_size' new objects to the stack.
*/
//inline CondLikelihoodStorage::CondLikelihoodStorage(unsigned cond_like_len, unsigned starting_size)
// :
//  cl_len(cond_like_len),
//  realloc_min(1)
//	{
//	assert(cl_len > 0);
//	if (starting_size > 0)
//		fillTo(starting_size);
//	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `cl_stack' is empty, calls CondLikelihoodStorage::fillTo to add `realloc_min' more objects to the stack. After
|	ensuring that the `cl_stack' is not empty, pops a CondLikelihoodShPtr off and returns it.
*/
inline CondLikelihoodShPtr CondLikelihoodStorage::getCondLikelihood()
	{
	assert(num_patterns > 0);
	assert(num_rates > 0);
	assert(num_states > 0);
	if (cl_stack.empty())
		fillTo(realloc_min);

	CondLikelihoodShPtr cl_ptr = cl_stack.top();
	cl_stack.pop();

#if defined(OBSOLETE_DEBUGGING_CODE)
	unsigned bytes = cl_ptr->getCLASize();
	unsigned expected_bytes = num_patterns*num_rates*num_states;
	if (bytes < expected_bytes)
		std::cerr << "***** bad, Bad, BAD! supplying CLA of length " << bytes << ", which is shorter than it needs to be (" << expected_bytes << ")" << std::endl;
#endif

	return cl_ptr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Pushes supplied CondLikelihoodShPtr `p' onto `cl_stack'.
*/
inline void CondLikelihoodStorage::putCondLikelihood(CondLikelihoodShPtr p)
	{
	assert(num_patterns > 0);
	assert(num_rates > 0);
	assert(num_states > 0);

#if defined(OBSOLETE_DEBUGGING_CODE)
	unsigned bytes = p->getCLASize();
	unsigned expected_bytes = num_patterns*num_rates*num_states;
	if (bytes < expected_bytes)
		std::cerr << "***** bad, Bad, BAD! storing CLA of length " << bytes << ", which is shorter than the current expected length (" << expected_bytes << ")" << std::endl;
#endif

	cl_stack.push(p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `realloc_min', which determines how many CondLikelihood objects are to be created when the 
|	`cl_stack' is discovered to be empty.
*/
inline void CondLikelihoodStorage::setReallocMin(unsigned sz)
	{
	assert(sz > 0);
	realloc_min = sz;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes all CondLikelihoodShPtr objects currently in the `cl_stack'.
*/
inline void CondLikelihoodStorage::clearStack()
	{
	while (!cl_stack.empty())
		cl_stack.pop();
	num_created = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data members `num_patterns', `num_rates' and `num_states', which determine the dimensions of all 
|	CondLikelihood objects stored. If the `cl_stack' is not currently empty and if the new conditional likelihood array
|	length is greater than the current length, all existing objects in `cl_stack' are deleted so that CLAs supplied to
|	the tree in the future will be at least the minimum length needed. Because it is critical that CLAs already checked
|	out to a tree be removed if the new length is longer than the old length (otherwise, CLAs that are too short will 
|	be used in the future), this function also checks to make sure there are no CLAs currently checked out. If there are
|	checked out CLAs, an XLikelihood exception is thrown. This function has no effect unless the supplied arguments
|	imply that newly-created CLAs will be longer than the existing ones. It is somewhat wastefull to leave in CLAs that
|	are longer than they need to be, but perhaps more wasteful to continually recall and delete all existing CLAs just
|	to ensure that their length is exactly correct.
*/
inline void CondLikelihoodStorage::setCondLikeDimensions(unsigned np, unsigned nr, unsigned ns)
	{
	assert(np > 0);
	assert(nr > 0);
	assert(ns > 0);
	unsigned newlen = np*nr*ns;
	unsigned oldlen = num_patterns*num_rates*num_states;

	if (newlen > oldlen)
		{
#if defined(OBSOLETE_DEBUGGING_CODE)
		unsigned checkouts = num_created - cl_stack.size();
		if (checkouts > 0)
			std::cerr << "### The number of among-site rates increased, but some shorter conditional likelihood arrays remain in the tree" << std::endl;
#endif
		clearStack();
		num_patterns = np;
		num_rates = nr;
		num_states = ns;
		}
	}

} //namespace phycas

#endif
