#if ! defined(COND_LIKELIHOOD_INL)
#define COND_LIKELIHOOD_INL

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
  underflowExponVec(npatterns),
  numEdgesSinceUnderflowProtection(UINT_MAX),
  claVec(npatterns*nrates*nstates)
	{
	cla = &claVec[0];
	uf = &underflowExponVec[0];
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
inline UnderflowType * CondLikelihood::getUF() const
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

} //namespace phycas

#endif
