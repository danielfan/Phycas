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

#include "mcmc_param.hpp"
#include "codon_model.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor and initializes its Codon pointer to NULL. Also sets the
|	`curr_value' data member to 1.0 and refreshes `curr_ln_prior' accordingly.
*/
OmegaParam::OmegaParam()
  : MCMCUpdater(), codon(NULL)
	{
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool OmegaParam::update()
	{
	//@POL should probably put next two lines in Model::update and have Model::update call a virtual Model::updateImpl
	// because it will be hard to remember to put these lines in every overloaded update function
	if (is_fixed)
		return false;
	slice_sampler->Sample();
    
	ChainManagerShPtr p = chain_mgr.lock();

    if (save_debug_info)
        {
        debug_info = str(boost::format("OmegaParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	return true;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class version adds the current omega parameter value to the data already stored in 
|	`fitting_sample'.
*/
void OmegaParam::educateWorkingPrior()
	{
	if (!isFixed())
		{
		PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
		double omega = getCurrValueFromModel();
		fitting_sample.push_back(omega);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use samples in `fitting_sample' to parameterize `ref_dist'. This function is called during s-cubed style
|	steppingstone sampling after the initial phase of sampling from the posterior so that the working prior can be
|	used for the remaining phases. Assumes `fitting_sample' has more than 1 element.  Assigns a GammaDistribution object
|	to `ref_dist'.
*/
void OmegaParam::finalizeWorkingPrior()
	{
	PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
	fitGammaWorkingPrior();
	}
}
