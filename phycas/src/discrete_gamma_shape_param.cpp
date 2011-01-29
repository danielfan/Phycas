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

#include "phycas/src/mcmc_param.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor. Also sets the `curr_value' data member to 0.5 (or 2.0
|	if `invert' is true. Sets `invert_shape' to the value of `invert'.
*/
DiscreteGammaShapeParam::DiscreteGammaShapeParam(bool invert)
  : MCMCUpdater(), shape_inverted(invert)
	{
    curr_value = (shape_inverted ? 2.0 : 0.5);
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Destructor.
*/
DiscreteGammaShapeParam::~DiscreteGammaShapeParam() 
	{
 	//std::cerr << "\n>>>>> DiscreteGammaShapeParam dying...likelihood use count = " << likelihood.use_count() << std::endl;
 	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool DiscreteGammaShapeParam::update()
	{
	if (is_fixed)
		{
		return false;
		}
		
	slice_sampler->Sample();
    
	ChainManagerShPtr p = chain_mgr.lock();
	
    if (save_debug_info)
        {
        debug_info = str(boost::format("DiscreteGammaShapeParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set the value of `gamma_shape' in `model' to `v'.
*/
void DiscreteGammaShapeParam::sendCurrValueToModel(double v)
	{
	if (shape_inverted)
		model->setShape(1.0/v);	// change the gamma shape parameter in the model
	else
		model->setShape(v);	// change the gamma shape parameter in the model
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return the value of `gamma_shape' in `model'.
*/
double DiscreteGammaShapeParam::getCurrValueFromModel() const
	{
	double a = model->getShape();
    return (shape_inverted ? (1.0/a) : a);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	DiscreteGammaShapeParam is a functor whose operator() returns a value proportional to the full-conditional posterior
|	probability density for a particular value of the gamma shape parameter. If the supplied gamma shape value `a' is
|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative
|	infinity).
*/
double DiscreteGammaShapeParam::operator()(
  double a)	/**< is a new value for the gamma shape parameter */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	if (a > 0.0)
		{
		sendCurrValueToModel(a);
		recalcPrior(); // base class function that recomputes curr_ln_prior for the value curr_value
		likelihood->recalcRelativeRates();	// must do this whenever model's shape parameter changes
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
		curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

		if (is_standard_heating)
			if (use_ref_dist)
				{
				PHYCAS_ASSERT(ref_dist);
				double curr_ln_ref_dist = ref_dist->GetLnPDF(a);
				return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_ref_dist;
				}
			else 
				return heating_power*(curr_ln_like + curr_ln_prior);
		else
			return heating_power*curr_ln_like + curr_ln_prior;
		}
    else
        return ln_zero;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class version adds the current gamma shape parameter value to the data already stored in 
|	`fitting_sample'.
*/
void DiscreteGammaShapeParam::educateWorkingPrior()
	{
	if (!isFixed())
		{
		PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
		double shape = getCurrValueFromModel();
		fitting_sample.push_back(shape);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use samples in `fitting_sample' to parameterize `ref_dist'. This function is called during s-cubed style
|	steppingstone sampling after the initial phase of sampling from the posterior so that the working prior can be
|	used for the remaining phases. Assumes `fitting_sample' has more than 1 element. Assigns a GammaDistribution object
|	to `ref_dist'.
*/
void DiscreteGammaShapeParam::finalizeWorkingPrior()
	{
	PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
	fitGammaWorkingPrior();
	//fitLognormalWorkingPrior();
	}
}
