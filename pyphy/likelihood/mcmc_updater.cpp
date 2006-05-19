#include "phycas/force_include.h"
#include "CipresCommlib/CipresDataMatrixHelper.h"
#include "pyphy/likelihood/likelihood_models.hpp"
#include "pyphy/likelihood/tree_likelihood.hpp"
//#include "phycas/rand/probability_distribution.hpp"
#include "pyphy/likelihood/mcmc_updater.hpp"
#include "pyphy/likelihood/mcmc_param.hpp"
//#include "pyphy/likelihood/mcmc_chain_manager.hpp"
#include "pyphy/phylogeny/basic_tree.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Calls `likelihood'->calcLnL(*`tree') to recalculate `curr_ln_like'. This function is important because the 
|	tempting getLnLike() member function only returns the value of `curr_ln_like' (it does not recalculate anything).
*/
double MCMCUpdater::recalcLike()
	{
	curr_ln_like = likelihood->calcLnL(*tree->GetFirstPreorder());

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
	//POL temporary!
	//if (dynamic_cast<EdgeLenMasterParam *>((MCMCUpdater *)this))
	//	std::cerr << "~~~~   " << curr_ln_prior << " (" << getName() << ")" << std::endl;

	return curr_ln_prior;
	}

}	// namespace phycas
