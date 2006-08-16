#if ! defined(COND_LIKELIHOOD_INL)
#define COND_LIKELIHOOD_INL

#include <numeric>
#include "boost/format.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	CondLikelihood constructor. Allocates `npatterns'*`nrates'*`nstates' elements to the conditional likelihood vector 
|	`claVec', allocates `npatterns' elements to the vector `underflowExpon', and sets `numEdgesSinceUnderflowProtection'
|	to UINT_MAX. Sets data member `cla' to point to the first element of the `claVec' vector and `uf' to point to the
|	first element of `underflowExponVec'. All three arguments shoudl be non-zero, and no error checking is done to 
|	ensure this because CondLikelihood objects are managed exclusively by CondLikelihoodStorage class, which ensures
|	that the dimensions are valid.
*/
inline CondLikelihood::CondLikelihood(
  unsigned npatterns,	/**< is the number of data petterns */
  unsigned nrates,		/**< is the number of among-site relative rate categories */
  unsigned nstates)		/**< is the number of states */
  :
  claVec(npatterns*nrates*nstates),
  uf_sum(0),
  underflowExponVec(npatterns),
  numEdgesSinceUnderflowProtection(UINT_MAX)
	{
	cla = &claVec[0];
	uf = &underflowExponVec[0];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the data member `numEdgesSinceUnderflowProtection'.
*/
inline unsigned	CondLikelihood::getUnderflowNumEdges() const
	{
	return numEdgesSinceUnderflowProtection;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of the data member `numEdgesSinceUnderflowProtection' to `n'.
*/
inline void	CondLikelihood::setUnderflowNumEdges(
  unsigned n)	/**< is the number of edges traversed since underflow correction was last applied */
	{
	numEdgesSinceUnderflowProtection = n;
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
|	Returns pointer to the underflow correction array stored by the vector data member `underflowExponVec'.
*/
inline UnderflowType * CondLikelihood::getUF()
	{
	return uf;	//PELIGROSO
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns pointer to the conditional likelihood array stored by the vector data member `underflowExponVec'.
*/
inline UnderflowType const * CondLikelihood::getUF() const
	{
	return uf;	//PELIGROSO
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns length of the vector data member `claVec'.
*/
inline unsigned CondLikelihood::getCLASize() const
	{
	return (unsigned)claVec.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns length of the vector data member `underflowExponVec'.
*/
inline unsigned CondLikelihood::getUFSize() const
	{
	return (unsigned)underflowExponVec.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets all elements in `underflowExponVec' to zero and sets `uf_sum' to 0 as well.
*/
inline void CondLikelihood::zeroUF()
	{
	uf_sum = 0;
	underflowExponVec.assign(underflowExponVec.size(), 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns sum of underflow corrections over all sites (i.e. value of data member `uf_sum').
*/
inline UnderflowType CondLikelihood::getUFSum() const
	{
	return uf_sum;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns reference to `uf_sum' data member.
*/
inline UnderflowType & CondLikelihood::getUFSumRef()
	{
	return uf_sum;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string containing a space-separated list of all elements in `underflowExponVec'.
*/
inline std::string CondLikelihood::debugShowUF() const
	{
	std::string s;
	for (std::vector<UnderflowType>::const_iterator it = underflowExponVec.begin(); it != underflowExponVec.end(); ++it)
		s += str(boost::format("%d ") % (*it));
	return s;
	}

} //namespace phycas

#endif
