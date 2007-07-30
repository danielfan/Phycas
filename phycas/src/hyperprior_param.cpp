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
|	The constructor of HyperPriorParam simply calls the base class (MCMCUpdater) constructor. Like other parameters,
|   sets `is_move' to false and `has_slice_sampler' to true. Because this is a hyperparameter, sets `is_hyper_param' to 
|   true. Also sets `is_master_param' to false and `curr_value' to 0.1.
*/
HyperPriorParam::HyperPriorParam()
  : MCMCUpdater()
	{
    //@POL this default constructor should never be used (but needs to be defined to avoid errors related to 
    // python/C++ container conversions)
    std::cerr << "*** fatal error: default constructor used for HyperPriorParam ***" << std::endl;
    std::exit(0);
	curr_value = 0.1;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor of HyperPriorParam simply calls the base class (MCMCUpdater) constructor. Like other parameters,
|   sets `is_move' to false and `has_slice_sampler' to true. Because this is a hyperparameter, sets `is_hyper_param' to 
|   true. Also sets `is_master_param' to false, `curr_value' to 0.1 and `edgelen_master_param' to the supplied edge
|   length master parameter shared pointer `p'.
*/
HyperPriorParam::HyperPriorParam(
  EdgeLenMasterParamShPtr p)    /**> is the edge length master parameter whose prior this parameter controls */
  : MCMCUpdater(), edgelen_master_param(p)
	{
	curr_value = 0.1;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
void HyperPriorParam::update()
	{
	if (is_fixed)
		return;
	if (slice_sampler)
		{
		slice_sampler->Sample();
		}

    if (save_debug_info)
        {
        if (slice_sampler)
            debug_info = str(boost::format("HyperPriorParam %f") % (slice_sampler->GetLastSampledXValue()));
        else
            debug_info = "HyperPriorParam (no slice sampler)";
        }
	}

}
