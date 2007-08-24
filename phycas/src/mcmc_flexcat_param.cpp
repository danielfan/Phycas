/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//#include "phycas/force_include.h"
#include "phycas/src/cipres/CipresDataMatrixHelper.h"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/mcmc_param.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/basic_tree.hpp"

namespace phycas
{

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set `curr_value' to the corresponding value in the `gamma_rates_unnorm' vector in 
|   `model'.
*/
void FlexRateParam::setCurrValueFromModel()
	{
    curr_value = model->getFlexRateUnnorm(which);
	}
#endif

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
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);
		}

#if POLPY_NEWWAY
    if (is_standard_heating)
        return heating_power*(curr_ln_like + curr_ln_prior);
    else
        return heating_power*curr_ln_like + curr_ln_prior;
#else
    return curr_ln_like + curr_ln_prior;
#endif
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set `curr_value' to the corresponding value in the `gamma_rate_probs' vector in 
|   `model'.
*/
void FlexProbParam::setCurrValueFromModel()
	{
    curr_value = model->getFlexProbUnnorm(which);
	}
#endif

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
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);
		}

#if POLPY_NEWWAY
    if (is_standard_heating)
        return heating_power*(curr_ln_like + curr_ln_prior);
    else
        return heating_power*curr_ln_like + curr_ln_prior;
#else
    return curr_ln_like + curr_ln_prior;
#endif
	}

}	// namespace phycas
