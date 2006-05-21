#if ! defined(COND_LIKELIHOOD_HPP)
#define COND_LIKELIHOOD_HPP

#include <vector>
#include <stack>
namespace phycas
{
typedef double LikeFltType;

class CondLikelihood
	{
	public:
		CondLikelihood(unsigned len)
			:underflowExpon(0),
			numEdgesSinceUnderflowProtection(UINT_MAX),
			claVec(len)
			{
			cla = &claVec[0];
			}

		/* enum CLADirection {CLA_DIR_UP, CLA_DIR_DOWN}; 
		
		CLADirection getDirection() const
			{
			return direction;
			}
		*/
		LikeFltType *	getCLA()
			{
			return cla;
			}
		LikeFltType *	getCLA() const
			{
			return cla;
			}
		unsigned	getCLASize() const
			{
			return claVec.size();
			}
	private:
			//CLADirection 	direction;
			LikeFltType *	cla;
			unsigned 		underflowExpon;
			unsigned 		numEdgesSinceUnderflowProtection;
			std::vector<LikeFltType> 	claVec;
	};

/*----------------------------------------------------------------------------------------------------------------------
|	A stack that stores pointers to CondLikelihood objects. CondLikelihoodStorage object can be asked for a pointer to a
|	CondLikelihood object when one needed for a likelihood calculation. If the stack is not empty, the pointer on top of
|	the stack is popped off and returned. If the stack is empty, several new CondLikelihood objects are created 
|	(`realloc_min' to be exact) on the heap and a pointer to one of them is returned. The `realloc_min' data member is 
|	by default 1 but can be modified to improved efficiency if such "stack faults" are expected to be common.
*/
class CondLikelihoodStorage
	{
	public:

										CondLikelihoodStorage(unsigned cond_like_len, unsigned starting_size = 1);

		CondLikelihood *				getCondLikelihood();
		void							fillTo(unsigned capacity);
		void							setReallocMin(unsigned sz);

	private:

		unsigned						cl_len;			/**< The length of a conditional likelihood array (needed for the CondLikelihood constructor) */
		unsigned						realloc_min;	/**< When a request is made and `cl_stack' is empty, `realloc_min' new objects are created and added to the stack */
		std::stack<CondLikelihood *>	cl_stack;		/**< The stack of CondLikelihood pointers */
	};

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
	assert(starting_size > 0);
	fillTo(starting_size);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `cl_stack' is empty, calls CondLikelihoodStorage::fillTo to add `realloc_min' more objects to the stack. After
|	ensuring that the `cl_stack' is not empty, pops a CondLikelihood pointer off and returns it.
*/
inline CondLikelihood * CondLikelihoodStorage::getCondLikelihood()
	{
	if (cl_stack.empty())
		fillTo(realloc_min);

	CondLikelihood * cl_ptr = cl_stack.top();
	cl_stack.pop();
	return cl_ptr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Ensures that the `cl_stack' contains at least `capacity' CondLikelihood pointers.
*/
inline void CondLikelihoodStorage::fillTo(unsigned capacity)
	{
	unsigned curr_sz = cl_stack.size();
	unsigned num_needed = (capacity > curr_sz ? capacity - curr_sz : 0);
	for (unsigned i = 0; i < num_needed; ++i)
		{
		cl_stack.push(new CondLikelihood(cl_len));
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `realloc_min', which determines how many CondLikelihood objects are to be added when the 
|	`cl_stack' is discovered to be empty.
*/
inline void CondLikelihoodStorage::setReallocMin(unsigned sz)
	{
	assert(sz > 0);
	realloc_min = sz;
	}

} //namespace phycas
#endif
