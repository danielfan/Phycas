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
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/mcmc_updater.hpp"
#include "phycas/src/mcmc_param.hpp"
#include "phycas/src/basic_tree.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|   Set `curr_value' data member to the supplied value `x'. Not ordinarily used, but useful for debugging.
*/
void MCMCUpdater::setCurrValue(
  double x) /**< is the value to which `curr_value' should be set */
	{
    curr_value = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing. Derived classes representing model parameters should override this method to 
|   set their `curr_value' data member to the corresponding value in `model'.
*/
void MCMCUpdater::setCurrValueFromModel()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls `likelihood'->calcLnL(`tree') to recalculate `curr_ln_like'. This function is important because the 
|	tempting getLnLike() member function only returns the value of `curr_ln_like' (it does not recalculate anything).
*/
double MCMCUpdater::recalcLike()
	{
	likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
	curr_ln_like = likelihood->calcLnL(tree);

	//std::cerr << "--> in MCMCUpdater::recalcLike(), name = " << getName();	//POL temp
	//std::cerr << ", curr_ln_like = " << curr_ln_like;	//POL temp
	//std::cerr << std::endl;	//POL temp

	return curr_ln_like;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the log prior just after this updater's update() member function was last called.
*/
double MCMCUpdater::getLnPrior() const
	{
	return curr_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the power used in heating (variable `heating_power') to the specified value `p'.
*/
void MCMCUpdater::setPower(
  double p)
	{
	heating_power = p;
    if (slice_sampler)
        {
        if (heating_power > 0.0)
            slice_sampler->UseDoublingMethod(false);
        else
            slice_sampler->UseDoublingMethod(true);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the power used in heating (variable `heating_power').
*/
double MCMCUpdater::getPower() const
	{
	return heating_power;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `is_standard_heating' to true.
*/
void MCMCUpdater::setStandardHeating()
	{
	is_standard_heating = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `is_standard_heating' to false.
*/
void MCMCUpdater::setLikelihoodHeating()
	{
	is_standard_heating = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the data member `is_standard_heating' is true.
*/
bool MCMCUpdater::isStandardHeating() const
	{
	return is_standard_heating;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the data member `is_standard_heating' is false.
*/
bool MCMCUpdater::isLikelihoodHeating() const
	{
	return !is_standard_heating;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if the data member `heating_power' is set to exactly 1.0.
*/
bool MCMCUpdater::isNoHeating() const
	{
	return (heating_power == 1.0);
	}
}	// namespace phycas
