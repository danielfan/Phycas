//#include "phycas/force_include.h"
#include "pyphy/src/cipres/CipresDataMatrixHelper.h"
#include "pyphy/src/likelihood_models.hpp"
#include "pyphy/src/tree_likelihood.hpp"
#include "pyphy/src/probability_distribution.hpp"
#include "pyphy/src/mcmc_param.hpp"
#include "pyphy/src/mcmc_chain_manager.hpp"
#include "pyphy/src/basic_tree.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	FlexRateParam is a functor whose operator() returns a value proportional to the full-conditional posterior
|	probability density for a particular value of an among-sites relative rate parameter in the FLEXCAT model. If the 
|	supplied relative rate `r' is out of bounds (i.e. < 0.0), the return value is `ln_zero' (closest we can come to a 
|	log posterior equal to negative infinity).
*/
double FlexRateParam::operator()(
  double r) /**< is a new value for the relative rate parameter */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	refreshLeftRightValues();

    if (r > left_value && r < right_value)
		{
		model->setFlexRateUnnorm(which, r);
		likelihood->recalcRelativeRates();	// must do this whenever one of the relative rates changes
        curr_value = r;
		recalcPrior();
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
        curr_ln_like = likelihood->calcLnL(tree);
		ChainManagerShPtr p = chain_mgr.lock();
		assert(p);
		p->setLastLnLike(curr_ln_like);
		}

	return curr_ln_like + curr_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	FlexProbParam is a functor whose operator() returns a value proportional to the full-conditional posterior
|	probability density for a particular value of an among-sites relative rate category probability parameter in the 
|	FLEXCAT model. If the supplied relative category probability `f' is out of bounds (i.e. < 0.0), the return value is 
|	`ln_zero' (closest we can come to a log posterior equal to negative infinity).
*/
double FlexProbParam::operator()(
  double f) /**< is a new value for the relative rate probability parameter */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

    if (f > 0.0)
		{
		model->setFlexProbUnnorm(which, f);
		likelihood->recalcRelativeRates();	// must do this whenever one of the rate probabilities changes
        curr_value = f;
		recalcPrior();
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
        curr_ln_like = likelihood->calcLnL(tree);
		ChainManagerShPtr p = chain_mgr.lock();
		assert(p);
		p->setLastLnLike(curr_ln_like);
		}

	return curr_ln_like + curr_ln_prior;
	}

}	// namespace phycas
