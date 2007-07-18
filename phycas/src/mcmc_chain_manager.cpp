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
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/basic_tree.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Refreshes `last_ln_like' by calling recalcLike() for the first updater in the `all_updaters' vector.
*/
void MCMCChainManager::refreshLastLnLike()
	{
	// Might need to add more checks here
	if (dirty)
		{
		throw XLikelihood("cannot call refreshLastLnLike() for chain manager before calling finalize()");
		}

	// Only need one updater to get the log-likelihood because all updaters must be able to do this calculation
	MCMCUpdaterShPtr u = *params_begin;
	last_ln_like = u->recalcLike();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Refreshes `last_ln_prior' by calling recalcPrior() for all parameters.
*/
void MCMCChainManager::refreshLastLnPrior()
	{
	// Might need to add more checks here
	if (dirty)
		{
		throw XLikelihood("cannot call refreshLastLnPrior() for chain manager before calling finalize()");
		}

	// Go through all parameters and get each of them to recalculate the prior probability density for their
	// current value (don't call getLnPrior() because this will return value of a data member but will not 
	// force a recalculation)
	last_ln_prior = 0.0;
	for (MCMCUpdaterIter p = params_begin; p != params_end; ++p)
		{
		MCMCUpdaterShPtr s = *p;
		double this_ln_prior = s->recalcPrior();
		last_ln_prior += this_ln_prior;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds edge length parameters (if `separate_edgelen_params' is true) and model-specific parameters to the 
|	`edge_len_params' and `model_params' data members, respectively. Also adds edge length hyperparameters to the 
|   `edge_len_hyperparams' data member if any hyperpriors have been specified. The model-specific parameters are added 
|   by calling the model's createParameters() member function. This function does not call the finalize() function, 
|   however, in case other parameters need to be added after this function is called.
*/
void MCMCChainManager::addMCMCUpdaters(
  ModelShPtr m,						/**< is the substitution model */
  TreeShPtr t,						/**< is the tree */
  TreeLikeShPtr like,				/**< is the likelihood calculator */
  LotShPtr r,						/**< is the pseudo-random number generator */
  unsigned max_units,				/**< is the maximum number of slice sampler units to use in each update */
  unsigned weight)					/**< is the weight to be used for all parameters added by this function */
	{
	if (!like)
		{
		throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no TreeLikelihood object defined");
		}
	if (!t)
		{
		throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no tree defined");
		}
	if (!r)
		{
		throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no pseudorandom number generator defined");
		}
	if (!m)
		{
		throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no model defined");
		}
	if (!m->getInternalEdgeLenPrior())
		{
		throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no internal edge length prior defined");
		}
	if (!m->getExternalEdgeLenPrior())
		{
		throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no external edge length prior defined");
		}

	MCMCUpdaterIter iter;
	MCMCUpdaterVect edgelens;
	MCMCUpdaterVect edgelen_hyperparams;
	MCMCUpdaterVect parameters;

	// Ask the model to create the edge length parameters, edge length hyperparameters, and 
	// its own model-specific parameters.
	m->createParameters(t, edgelens, edgelen_hyperparams, parameters);

	// Add the edge length parameters (might be master parameters for internal and external edges or, 
    // if separate_edgelen_params is true, one edge length parameter for every edge in the tree)
	for (iter = edgelens.begin(); iter != edgelens.end(); ++iter)
		{
		MCMCUpdaterShPtr p = (*iter);
		p->setWeight(weight);
		p->setMaxUnits(max_units);
		p->setModel(m);
		p->setTreeLikelihood(like);
		p->setLot(r);
		addEdgeLenParam(p);
		}

	// Add the edge length hyperparameters (if any were created)
	for (iter = edgelen_hyperparams.begin(); iter != edgelen_hyperparams.end(); ++iter)
		{
		MCMCUpdaterShPtr p = (*iter);
		p->setWeight(weight);
		p->setMaxUnits(max_units);
		p->setModel(m);
		p->setTreeLikelihood(like);
		p->setLot(r);
		addEdgeLenHyperparam(p);
		}

	for (iter = parameters.begin(); iter != parameters.end(); ++iter)
		{
		MCMCUpdaterShPtr p = (*iter);
		p->setWeight(weight);
		p->setMaxUnits(max_units);
		p->setModel(m);
		p->setTreeLikelihood(like);
		p->setLot(r);
		addModelParam(p);
		}
	}

}	// namespace phycas
