#if ! defined(COND_LIKELIHOOD_HPP)
#define COND_LIKELIHOOD_HPP

#include <vector>
#include <stack>
namespace phycas
{
typedef double LikeFltType;
typedef unsigned UnderflowType;

/*----------------------------------------------------------------------------------------------------------------------
|	Manages a conditional likelihood array for one end of an edge.
*/
class CondLikelihood
	{
	public:

									CondLikelihood(unsigned npatterns, unsigned nrates, unsigned nstates);

		LikeFltType *				getCLA();
		LikeFltType *				getCLA() const;
		unsigned					getCLASize() const;

		UnderflowType *				getUF();
		UnderflowType const *		getUF() const;
		unsigned					getUFSize() const;

		void						zeroUF();
		std::string					debugShowUF() const;
		UnderflowType				getUFSum() const;
		UnderflowType &				getUFSumRef();

		unsigned					getUnderflowNumEdges() const;
		void						setUnderflowNumEdges(unsigned n);

	private:

		LikeFltType *				cla;								/**< Pointer to conditional likelihood array stored by `claVec'  */
		std::vector<LikeFltType>	claVec;								/**< Each element contains the likelihood conditional on a particular state, rate and pattern */

		UnderflowType				uf_sum;								/**< Stores sum of underflow corrections over all sites */
		UnderflowType *				uf;									/**< Pointer to the underflow correction array stored by `underflowExponVec'  */
		std::vector<UnderflowType>	underflowExponVec;					/**< Stores log of the underflow correction factor for each pattern */

		unsigned 					numEdgesSinceUnderflowProtection;	/**< The number of edges traversed since the underflow protection factor was last updated */
	};

typedef boost::shared_ptr<CondLikelihood> CondLikelihoodShPtr;
typedef boost::shared_ptr<const CondLikelihood> ConstCondLikelihoodShPtr; 

} //namespace phycas

#include "pyphy/src/cond_likelihood.inl"

#endif