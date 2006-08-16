#if ! defined(MCMC_MOVE_INL)
#define MCMC_MOVE_INL

#include <limits>		// for std::numeric_limits

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor for the base class of all move classes used for MCMC. MCMCMove-derived objects manage Metropolis-
|	Hastings proposals in a Bayesian analysis.
*/
inline MCMCMove::MCMCMove() : curr_ln_prior(0.0), curr_ln_like(0.0)
	{
	//ln_zero = std::log(std::numeric_limits<double>::denorm_min());
	ln_zero = -DBL_MAX;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a Metropolis-Hastings update, which involves calling proposeNewState, computing the posterior density of
|	the new state, and deciding whether to accept or reject the proposed move.
*/
inline void MCMCMove::update()
	{
	// TODO: let MCMCMove objects be stored in the MCMCChainManager
	// TODO: give each MCMCMove or MCMCParam stored in the MCMCChainManager an unsigned weight that equals the number
	//       of times to update that move or parameter 
	// TODO: always store last calculated posterior density in MCMCChainManager for easy reference
	// MCMCChainManager will eventually encapsulate an individual chain (of perhaps several coupled chains)
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The MCMCChainManager::finalize function calls this function for each MCMCMove in knows about. This provides a way 
|	for each move to call the MCMCChainManager when it needs the joint prior over all parameters.
*/
inline void MCMCMove::setChainManager(
  ChainManagerWkPtr p)		/**< is a pointer to the MCMCChainManager containing this parameter */
	{
	chain_mgr = p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the log prior just after this move's accept() or revert() member function was last called.
*/
inline double MCMCMove::getLnPrior() const
	{
	//std::cerr << "--> MCMCMove::getLnPrior(), name = " << getName() << std::endl;	//POL temp

	return curr_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the log likelihood just after this move's accept() or revert() member function was last called.
*/
inline double MCMCMove::getLnLike() const
	{
	return curr_ln_like;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a const reference to the `name' data member.
*/
inline const std::string & MCMCMove::getName() const
	{
	return name;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `name' data member of this move to the supplied string.
*/
inline void MCMCMove::setName(const std::string & s)
	{
	name = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `rng' data member to the supplied Lot shared pointer.
*/
inline void MCMCMove::setLot(LotShPtr r)
	{
	rng = r;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `tree' data member to the supplied Tree shared pointer.
*/
inline void MCMCMove::setTree(TreeShPtr p)
	{
	tree = p;
	tree_manipulator.setTree(p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `likelihood' data member to the supplied TreeLikelihood shared pointer.
*/
inline void MCMCMove::setTreeLikelihood(TreeLikeShPtr L)
	{
	likelihood = L;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `model' data member to the supplied Model shared pointer.
*/
inline void MCMCMove::setModel(ModelShPtr m)
	{
	model = m;
	}

} // namespace phycas
#endif
