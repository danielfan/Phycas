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
|	The StateFreqParam constructor requires the caller to specify a value for `which'. `which' is 0, 1, 2, or 3 for 
|	nucleotide models and  determines the particular state frequency (e.g. A, C, G, or T, respectively) this object 
|	represents. It calls the base class (MCMCUpdater) constructor. Also sets the `curr_value' data member to 1.0 and 
|	refreshes `curr_ln_prior' accordingly.
*/
StateFreqParam::StateFreqParam(
  unsigned w)		/**< The 0-based index of the base frequency being managed by this object (0=A, 1=C, 2=G and 3=T) */
  : MCMCUpdater(),  which(w)
	{
	curr_value = 1.0;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool StateFreqParam::update()
	{
	if (is_fixed)
		return false;
	slice_sampler->Sample();

    if (save_debug_info)
        {
        debug_info = str(boost::format("StateFreqParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	return true;
	}
}
