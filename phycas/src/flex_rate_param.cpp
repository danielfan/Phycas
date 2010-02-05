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
|	The FlexRateParam constructor requires a reference to an unsigned variable that holds the number of spacers (fake
|	rates that serve to determine the strength of the prior for relative rate parameters). It also requires a reference
|	to a double variable that holds the maximum possible value for a relative rate. Finally, a reference to the vector
|	of relative rates is required because in order to update a particular rate parameter, one needs to know the rates 
|	on either side of it.
*/
FlexRateParam::FlexRateParam(
  unsigned & s,					/**< is the number of spacers between each rate in the order statistics prior */
  double & ub,					/**< is the maximum possible value for a rate */
  std::vector<double> & rr)		/**< is the relative rate parameter vector (this parameter needs to know the value on either side when updating) */
  : MCMCUpdater(),  nspacers(s), upper_bound(ub), 
  flex_rates(rr),
  left_value(0.0), right_value(1.0), which(0)
	{
	curr_value = flex_rates[0];
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The MCMCChainManager::finalize function calls this function for each MCMCUpdater it knows about. This provides a way 
|	for each updater to call back the MCMCChainManager when it needs the joint prior over all parameters.
*/
void FlexRateParam::setChainManager(
 ChainManagerWkPtr p)		/**< is a pointer to the MCMCChainManager containing this updater */
	{
	MCMCUpdater::setChainManager(p);
	//refreshLeftRightValues();
	//double starting_x = (left_value + right_value)/2.0;
	//slice_sampler->SetXValue(starting_x);
	//slice_sampler->SetXValue(curr_value);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Determines the quantities `curr_value', `left_value' and `right_value'. The smallest this parameter can become is 
|	`left_value', which should equal `flex_rates'[`which' - 1] (unless `which' is 0, in which case `left_value' is 0.0).
|	The largest this parameter can become is `right_value', which should equal `flex_rates'[`which' + 1] (unless `which'
|	is `num_rates' - 1, in which case `right_value' is `upper_bound'). The `curr_value' of this parameter is equal to
|	`flex_rates'[`which'].
*/
void FlexRateParam::refreshLeftRightValues()
	{
	PHYCAS_ASSERT(flex_rates.size() > 0);
	unsigned last = (unsigned)flex_rates.size() - 1;
	if (which == 0 && which == last)
		{
		// only one rate (i.e. no rate heterogeneity)
		left_value  = 0.0;
		right_value = upper_bound;
		}
	else if (which == 0)
		{
		left_value = 0.0;
		right_value = flex_rates[which + 1];
		}
	else if (which == last)
		{
		left_value = flex_rates[which - 1];
		right_value = upper_bound;
		}
	else
		{
		left_value = flex_rates[which - 1];
		right_value = flex_rates[which + 1];
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Order statistics prior" (user gets no choice of prior in this case).
*/
std::string FlexRateParam::getPriorDescr() const
	{
	return std::string("Order statistics prior");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the order statistics prior used for relative rates in the FLEX model, and sets `curr_ln_prior'. The order 
|	statistics prior is largest when `curr_value' is between `left_value' and `right_value' and becomes tiny if 
|	`curr_value' gets close to either `left_value' or `right_value'. Between each pair of rates lies one or more 
|   imaginary spacers. The more spacers lying between each adjacent pair of rates, the harder it is for rates to get 
|   close together. The number of spacers between each adjacent pair is specified by the data member `nspacers', which 
|   is set in the constructor.
*/
double FlexRateParam::recalcPrior()
	{
	refreshLeftRightValues();
#if POLPY_NEWWAY
    double curr_value = flex_rates[which];
#else //old way
    PHYCAS_ASSERT(curr_value == flex_rates[which]);
#endif
	curr_ln_prior = (double)nspacers*(std::log(curr_value - left_value) + std::log(right_value - curr_value));
	return curr_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool FlexRateParam::update()
	{
	if (is_fixed)
		return false;

    // Ok, it looks like the following cannot work, but it does! The key is the assignment to `which'
    // and the fact that SliceSampler::Sample() calls FlexRateParam::operator() repeatedly. Inside
    // the operator() function, `which' is used to determine which of the relative rates is 
    // modified. See FlexRateParam::operator() in mcmc_flexcat_param.cpp for details.
	unsigned nrates = flex_rates.size();
	for (unsigned i = 0; i < nrates; ++i)
		{
		which = i;
		curr_value = flex_rates[which];
		slice_sampler->SetXValue(curr_value);
		slice_sampler->Sample();
		}

    if (save_debug_info)
        {
        debug_info = "FlexRateParam ";
    	for (unsigned i = 0; i < nrates; ++i)
            debug_info += str(boost::format("\n  %d %f") % i % flex_rates[i]);
        }
	return true;
	}
}
