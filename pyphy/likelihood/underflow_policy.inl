#if !defined(UNDERFLOW_POLICY_INL)
#define UNDERFLOW_POLICY_INL

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of the data member `underflow_num_edges' to `nedges'. This is the number of edges to accumulate before 
|	correcting for underflow. Assumes `nedges' is greater than zero. Often several hundred taxa are required before
|	underflow becomes a problem, so a reasonable value for `underflow_num_edges' is 25.
*/
inline void SimpleUnderflowPolicy::setTriggerSensitivity(
  unsigned nedges)		/**< is the number of edges to traverse before correcting for underflow */
	{
	assert(nedges > 0);
	underflow_num_edges = nedges;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of data member `underflow_max_value' to `maxval'. The largest conditional likelihood for every pattern 
|	will be adjusted to a value close to this. Suppose the largest conditional likelihood for pattern i was x. A factor 
|	exp(f) is found such that x*exp(f) = `underflow_max_value'. The fractional component of f is then removed, yielding
|	an unsigned int k, which is what is actually stored. The largest conditional likelihood for pattern i after the
|	correction is thus x*exp(k), which will be less than or equal to `underflow_max_value'. Assumes `maxval' is greater
|	than zero. A reasonable value for `underflow_max_value' is 10000.
*/
inline void SimpleUnderflowPolicy::setCorrectToValue(
  double maxval)	/**< is the target value to which the largest conditional likelihood for any given pattern will be scaled */
	{
	assert(maxval > 0.0);
	underflow_max_value = maxval;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `num_rates', `num_patterns', and `num_states' data members to `nr', `np' and `ns', respectively. Assumes
|	that all three values are greater than zero.
*/
inline void SimpleUnderflowPolicy::setDimensions(
  unsigned np,	/**< is the number of patterns in the data */
  unsigned nr, 	/**< is the number of relative rates used in modeling among-site rate variation */
  unsigned ns)	/**< is the number of states */
	{
	assert(np > 0);
	assert(nr > 0);
	assert(ns > 0);
	num_patterns = np;
	num_rates = nr;
	num_states = ns;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Handles case of a node subtending two tips. In this case, the number of edges traversed is just two, and thus no
|	corrective action is taken other than to set the number of underflow edges traversed to 2. The zeroUF member
|	function of the supplied `cond_like' object is called, however, to ensure that the cumulative correction factor is 
|	zero (this `cond_like' might have just been brought in from storage with some correction factor already in place).
*/
inline void SimpleUnderflowPolicy::twoTips(
  CondLikelihood & cond_like) /**< is the conditional likelihood array to correct */
  const
	{
	cond_like.setUnderflowNumEdges(2);
	cond_like.zeroUF();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Nothing needs to be done in the simple case, so this method is a no-op.
*/
inline void SimpleUnderflowPolicy::correctSiteLike(
  double & site_like,						/**< is the log-site-likelihood to correct */
  unsigned pat,								/**< is the index of the pattern representing the site to correct */
  ConstCondLikelihoodShPtr condlike_shptr)	/**< is the conditional likelihood array of the likelihood root node */
  const
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains a pointer to the underflow array for `cond_like', then subtracts the value stored in that underflow array 
|	for pattern `pat' from the site likelihood `site_like'.
*/
inline void PatternSpecificUnderflowPolicy::correctSiteLike(
  double & site_like,						/**< is the log-site-likelihood to correct */
  unsigned pat,								/**< is the index of the pattern representing the site to correct */
  ConstCondLikelihoodShPtr condlike_shptr)	/**< is the conditional likelihood array of the likelihood root node */
  const
	{
	assert(condlike_shptr);
	UnderflowType const * uf = condlike_shptr->getUF();
	assert(uf != NULL);
	site_like -= (double)uf[pat];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains sum of underflow correction for all patterns in `cond_like', then subtracts that value from the 
|	log-likelihood `ln_like'.
*/
inline void SimpleUnderflowPolicy::correctLnLike(
  double & ln_like,							/**< is the log-likelihood to correct for underflow */
  ConstCondLikelihoodShPtr condlike_shptr)	/**< is the conditional likelihood array of the likelihood root node */
  const
	{
	assert(condlike_shptr);
	UnderflowType ufsum = condlike_shptr->getUFSum();
	ln_like -= (double)ufsum;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Nothing needs to be done to the overall log-likelihood because the correction takes place at the site log-likelihood
|	stage, so this method is a no-op.
*/
inline void PatternSpecificUnderflowPolicy::correctLnLike(
  double & ln_like,							/**< is the log-likelihood to correct for underflow */
  ConstCondLikelihoodShPtr condlike_shptr)	/**< is the conditional likelihood array of the likelihood root node */
  const
	{
	}

} // namespace phycas

#endif
