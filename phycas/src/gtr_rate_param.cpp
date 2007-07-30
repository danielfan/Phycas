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
|	Constructor calls the base class (MCMCUpdater) constructor and initializes its GTR pointer to NULL. Also sets the
|	`curr_value' data member to 1.0 and refreshes `curr_ln_prior' accordingly. Assumes `w' is greater than or equal to
|	zero and less than 6.
*/
GTRRateParam::GTRRateParam(
  unsigned w)		/**< The 0-based index of the relative rate being managed by this object (0=AC, 1=AG, 2=AT, 3=CG, 4=CT and 5=GT) */
  : MCMCUpdater(), gtr(NULL), which(w)
	{
	PHYCAS_ASSERT(w >= 0);
	PHYCAS_ASSERT(w < 6);
	curr_value = 1.0;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
void GTRRateParam::update()
	{
	//@POL should probably put next two lines in Model::update and have Model::update call a virtual Model::updateImpl
	// because it will be hard to remember to put these lines in every overloaded update function
	if (is_fixed)
		return;
	slice_sampler->Sample();
    
    if (save_debug_info)
        {
        debug_info = str(boost::format("GTRRateParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	}

}
