#if ! defined(COND_LIKELIHOOD_STORAGE_INL)
#define COND_LIKELIHOOD_STORAGE_INL

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
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data members `num_patterns', `num_rates' and `num_states', which determine the dimensions of all 
|	CondLikelihood objects stored. If the `cl_stack' is not currently empty and if any one of `np', `nr' or `ns' differ
|	from their corresponding CondLikelihoodStorage data members, all existing objects in `cl_stack' are deleted.
*/
inline void CondLikelihoodStorage::setCondLikeDimensions(unsigned np, unsigned nr, unsigned ns)
	{
	assert(np > 0);
	assert(nr > 0);
	assert(ns > 0);
	bool changed = ((np != num_patterns) || (nr != num_rates) || (ns != num_states));
	if (changed)
		clearStack();
	num_patterns = np;
	num_rates = nr;
	num_states = ns;
	}

} //namespace phycas

#endif
