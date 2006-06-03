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

									CondLikelihood(unsigned len);

		LikeFltType *				getCLA();
		LikeFltType *				getCLA() const;
		unsigned					getCLASize() const;

	private:

		LikeFltType *				cla;								/**< Pointer to conditional likelihood array stored by `claVec'  */
		unsigned 					underflowExpon;						/**< Log of the underflow correction factor */
		unsigned 					numEdgesSinceUnderflowProtection;	/**< The number of edges traversed since the underflow protection factor was last updated */
		std::vector<LikeFltType>	claVec;								/**< Each element contains the likelihood conditional on a particular state, rate and pattern */
	};

typedef boost::shared_ptr<CondLikelihood> CondLikelihoodShPtr;
typedef boost::shared_ptr<const CondLikelihood> ConstCondLikelihoodShPtr; 

/*----------------------------------------------------------------------------------------------------------------------
|	A stack that stores CondLikelihood shared pointers. CondLikelihoodStorage object can be asked for a pointer to a
|	CondLikelihoodShPtr object when one needed for a likelihood calculation. If the stack is not empty, the pointer on 
|	top of the stack is popped off and returned. If the stack is empty, several new CondLikelihood objects are created 
|	(`realloc_min' to be exact) on the heap and a pointer to one of them is returned. The `realloc_min' data member is 
|	by default 1 but can be modified to improved efficiency if such "stack faults" are expected to be common.
*/
class CondLikelihoodStorage
	{
	public:

										CondLikelihoodStorage();
										CondLikelihoodStorage(unsigned cond_like_len, unsigned starting_size = 1);

										~CondLikelihoodStorage();

		CondLikelihoodShPtr				getCondLikelihood();
		void							putCondLikelihood(CondLikelihoodShPtr p);
		void							fillTo(unsigned capacity);

		void							setCondLikeLength(unsigned sz);
		void							setReallocMin(unsigned sz);
		void							clearStack();

	private:

		unsigned						cl_len;			/**< The length of a conditional likelihood array (needed for the CondLikelihood constructor) */
		unsigned						realloc_min;	/**< When a request is made and `cl_stack' is empty, `realloc_min' new objects are created and added to the stack */
		std::stack<CondLikelihoodShPtr>	cl_stack;		/**< The stack of CondLikelihoodShPtr */
	};

} //namespace phycas

#include "pyphy/likelihood/cond_likelihood.inl"

#endif
