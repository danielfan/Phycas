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
    
    if (save_debug_info)
        {
        debug_info = str(boost::format("DiscreteGammaShapeParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class version adds the current gamma shape parameter value to the data already stored in 
|	`fitting_sample'.
*/
void DiscreteGammaShapeParam::educateWorkingPrior()
	{
	PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
	double shape = getCurrValueFromModel();
	fitting_sample.push_back(shape);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use samples in `fitting_sample' to parameterize `working_prior'. This function is called during s-cubed style
|	steppingstone sampling after the initial phase of sampling from the posterior so that the working prior can be
|	used for the remaining phases. Assumes `fitting_sample' has more than 1 element. Assigns a GammaDistribution object
|	to `working_prior'.
*/
void DiscreteGammaShapeParam::finalizeWorkingPrior()
	{
	PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
	fitGammaWorkingPrior();
	}
}
