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
|	The FlexProbParam constructor requires a reference to the vector of category probabilities so that it can deposit
|	a changed value in the appropriate place after an update.
*/
FlexProbParam::FlexProbParam(
  std::vector<double> & rp)		/**< is the category probability parameter vector */
  : MCMCUpdater(), which(0), rate_probs(rp)
	{
	curr_value = rate_probs[0];
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

void FlexProbParam::setLot(LotShPtr r)
	{
	rng = r;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member to update each of the probabilities in the 
|   vector `rate_probs'.
*/
void FlexProbParam::update()
	{
	if (is_fixed)
		return;

    // Ok, it looks like the following cannot work, but it does! The key is the assignment to `which'
    // and the fact that SliceSampler::Sample() calls FlexProbParam::operator() repeatedly. Inside
    // the operator() function, `which' is used to determine which of the rate probabilities is 
    // modified. See FlexProbParam::operator() in mcmc_flexcat_param.cpp for details.
	unsigned nrates = rate_probs.size();
	for (unsigned i = 0; i < nrates; ++i)
		{
		which = i;
		curr_value = rate_probs[which];
		slice_sampler->SetXValue(curr_value);
		slice_sampler->Sample();
		}

    if (save_debug_info)
        {
        debug_info = "FlexProbParam ";
    	for (unsigned i = 0; i < nrates; ++i)
            debug_info += str(boost::format("\n  %d %f") % i % rate_probs[i]);
        }
	}

}
