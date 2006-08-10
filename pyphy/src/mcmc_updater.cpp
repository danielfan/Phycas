//#include "phycas/force_include.h"
#include "pyphy/src/cipres/CipresDataMatrixHelper.h"
#include "pyphy/src/likelihood_models.hpp"
#include "pyphy/src/tree_likelihood.hpp"
#include "pyphy/src/mcmc_updater.hpp"
#include "pyphy/src/mcmc_param.hpp"
#include "pyphy/src/basic_tree.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Calls `likelihood'->calcLnL(`tree') to recalculate `curr_ln_like'. This function is important because the 
|	tempting getLnLike() member function only returns the value of `curr_ln_like' (it does not recalculate anything).
*/
double MCMCUpdater::recalcLike()
	{
	likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
	curr_ln_like = likelihood->calcLnL(tree);

	//std::cerr << "--> in MCMCUpdater::recalcLike(), name = " << getName();	//POL temp
	//std::cerr << ", curr_ln_like = " << curr_ln_like;	//POL temp
	//std::cerr << std::endl;	//POL temp

	return curr_ln_like;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the log prior just after this updater's update() member function was last called.
*/
double MCMCUpdater::getLnPrior() const
	{
	return curr_ln_prior;
	}

}	// namespace phycas