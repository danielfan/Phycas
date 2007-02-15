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
#include "phycas/src/mcmc_updater.hpp"
#include "phycas/src/mcmc_param.hpp"
#include "phycas/src/basic_tree.hpp"

namespace phycas
{

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

}	// namespace phycas
