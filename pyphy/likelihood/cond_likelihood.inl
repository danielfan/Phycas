#if ! defined(COND_LIKELIHOOD_INL)
#define COND_LIKELIHOOD_INL

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	CondLikelihood constructor. Allocates `len' elements to the conditional likelihood vector `claVec', sets 
|	`numEdgesSinceUnderflowProtection' to UINT_MAX, and `underflowExpon' to 0.
*/
inline CondLikelihood::CondLikelihood(
  unsigned len)	/**< is the number of elements to allocate for `claVec' */
  :
  underflowExpon(0),
  numEdgesSinceUnderflowProtection(UINT_MAX),
  claVec(len)
	{
	cla = &claVec[0];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns pointer to the conditional likelihood array stored by the vector data member `claVec'.
*/
inline LikeFltType * CondLikelihood::getCLA()
	{
	return cla;	//PELIGROSO
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns pointer to the conditional likelihood array stored by the vector data member `claVec'.
*/
inline LikeFltType * CondLikelihood::getCLA() const
	{
	return cla;	//PELIGROSO
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns length of the vector data member `claVec'.
*/
inline unsigned CondLikelihood::getCLASize() const
	{
	return (unsigned)claVec.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `realloc_min' to 1 and `cl_len' to 0. 
*/
inline CondLikelihoodStorage::CondLikelihoodStorage()
  : 
  cl_len(0), 
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
inline CondLikelihoodStorage::CondLikelihoodStorage(unsigned cond_like_len, unsigned starting_size)
  :
  cl_len(cond_like_len),
  realloc_min(1)
	{
	assert(cl_len > 0);
	if (starting_size > 0)
		fillTo(starting_size);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `cl_stack' is empty, calls CondLikelihoodStorage::fillTo to add `realloc_min' more objects to the stack. After
|	ensuring that the `cl_stack' is not empty, pops a CondLikelihoodShPtr off and returns it.
*/
inline CondLikelihoodShPtr CondLikelihoodStorage::getCondLikelihood()
	{
	assert(cl_len > 0);
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
	assert(cl_len > 0);
	cl_stack.push(p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Ensures that the `cl_stack' contains at least `capacity' CondLikelihoodShPtr objects.
*/
inline void CondLikelihoodStorage::fillTo(unsigned capacity)
	{
	assert(cl_len > 0);
	unsigned curr_sz = (unsigned)cl_stack.size();
	unsigned num_needed = (capacity > curr_sz ? capacity - curr_sz : 0);
	for (unsigned i = 0; i < num_needed; ++i)
		{
		cl_stack.push(CondLikelihoodShPtr(new CondLikelihood(cl_len)));
		}
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
|	Sets the data member `cl_len', which determines the length of all CondLikelihood objects stored. If the	`cl_stack'
|	is not currently empty and if `sz' is not equal to the current value of `cl_len', deletes all existing objects in
|	`cl_stack'.
*/
inline void CondLikelihoodStorage::setCondLikeLength(unsigned sz)
	{
	assert(sz > 0);
	if (sz != cl_len)
		clearStack();
	cl_len = sz;
	}

} //namespace phycas

#endif
