/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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
#include "phypy/src/likelihood_models.hpp"
#include "phypy/src/mcmc_chain_manager.hpp"
#include "phypy/src/xlikelihood.hpp"
#include "phypy/src/basic_tree.hpp"

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
|	`edge_len_params' and `model_params' data members, respectively. Also sets `edge_len_hyperparam' if a prior for one
|	has been specified in the model. The model-specific parameters are added by calling the model's createParameters() 
|	member function. Does not call the finalize() function, however, in case other parameters need to be added after 
|	this function is called.
*/
void MCMCChainManager::addMCMCUpdaters(
  ModelShPtr m,						/**< is the substitution model */
  TreeShPtr t,						/**< is the tree */
  TreeLikeShPtr like,				/**< is the likelihood calculator */
  LotShPtr r,						/**< is the pseudo-random number generator */
  bool separate_edgelen_params,		/**< is true if one EdgeLenParam should be added for every edge and false if a single EdgeLenMasterParam should be added  */
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
	if (!m->getEdgeLenPrior())
		{
		throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no edge length prior defined");
		}

	MCMCUpdaterIter iter;
	MCMCUpdaterVect edgelens;
	MCMCUpdaterShPtr edgelen_hyperparam;
	MCMCUpdaterVect parameters;

	// Ask the model to create the edge length parameters, edge length hyperparameter, and 
	// its own model-specific parameters. The model is responsible for calling setName, setTree
	// and setPrior for each parameter it creates.
	m->createParameters(t, edgelens, edgelen_hyperparam, parameters, separate_edgelen_params);

	// Add the edge length parameters (might be a single edge length master parameter 
	// or, if separate_edgelen_params is true, one edge length parameter for every edge in the tree)
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

	// Add the edge length hyperparameter (if one was created)
	if (edgelen_hyperparam)
		{
		edgelen_hyperparam->setWeight(weight);
		edgelen_hyperparam->setMaxUnits(max_units);
		edgelen_hyperparam->setModel(m);
		edgelen_hyperparam->setTreeLikelihood(like);
		edgelen_hyperparam->setLot(r);
		addEdgeLenHyperparam(edgelen_hyperparam);
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
