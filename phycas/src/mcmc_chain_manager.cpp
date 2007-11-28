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
    // one edge length master parameter that governs all edges in the tree)
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

////////////////////////////////////////// inserted contents of mcmc_chain_manager.inl here /////////////////////////////

/*----------------------------------------------------------------------------------------------------------------------
|	The MCMCChainManager constructor does nothing.
*/
MCMCChainManager::MCMCChainManager() 
  : last_ln_like(0.0), last_ln_prior(0.0), dirty(true)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `all_updaters' data member as a const MCMCUpdaterVect reference. Used so that Python programs can iterate 
|	through all updaters using a construct such as this: for x in chainMgr.getAllUpdaters(): ...
*/
const MCMCUpdaterVect & MCMCChainManager::getAllUpdaters() const 
	{
	return all_updaters;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `moves' data member as a const MCMCUpdaterVect reference. Used so that Python programs can iterate 
|	through the stored moves using a construct such as this: for x in chainMgr.getMoves(): ...
*/
const MCMCUpdaterVect & MCMCChainManager::getMoves() const 
	{
	return moves;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `model_params' data member as a const MCMCUpdaterVect reference. Used so that Python programs can iterate 
|	through the model-specific parameters using a construct such as this: for x in chainMgr.getModelParams(): ...
*/
const MCMCUpdaterVect & MCMCChainManager::getModelParams() const 
	{
	return model_params;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `edge_len_params' data member as a const MCMCUpdaterVect reference. Used so that Python programs can iterate 
|	through the edge length parameters using a construct such as this: for x in chainMgr.getEdgeLenParams(): ...
*/
const MCMCUpdaterVect & MCMCChainManager::getEdgeLenParams() const 
	{
	return edge_len_params;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `edge_len_hyperparams' data member, which is a shared pointer to a HyperPriorParam.
|	Returns the `edge_len_hyperparams' data member as a const MCMCUpdaterVect reference. Used so that Python programs 
|   can iterate through the edge length hyperparameters using a construct such as this: 
|   for x in chainMgr.getEdgeLenHyperparams(): ...
*/
MCMCUpdaterVect MCMCChainManager::getEdgeLenHyperparams() const
	{
	return edge_len_hyperparams;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds a pointer to an MCMCUpdater-derived object to the `moves' vector. The `moves' vector is a vector of shared 
|	pointers to ensure that each move object is eventually destroyed.
*/
void	MCMCChainManager::addMove(MCMCUpdaterShPtr p) 
	{
	moves.push_back(p);
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds a pointer to an MCMCUpdater-derived object to the `model_params' vector. The `model_params' vector is a 
|	vector of shared pointers to ensure that each parameter object is eventually destroyed.
*/
void	MCMCChainManager::addModelParam(MCMCUpdaterShPtr p) 
	{
	model_params.push_back(p);
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds an EdgeLenParam shared pointer to the `edge_len_params' vector.
*/
void	MCMCChainManager::addEdgeLenParam(MCMCUpdaterShPtr p) 
	{
	edge_len_params.push_back(p);
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds a HyperPriorParam shared pointer to the `edge_len_hyperparams' vector.
*/
void	MCMCChainManager::addEdgeLenHyperparam(MCMCUpdaterShPtr p) 
	{
	edge_len_hyperparams.push_back(p);
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current log likelihood. Each call to an update() member function of one of the updaters refreshes the
|	`last_ln_like' data member.
*/
double MCMCChainManager::getLastLnLike() const 
	{
	return last_ln_like;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Used by an MCMCUpdater-derived object to set the current log prior just before it returns from its update()
|	function. Ordinarily, the calcJointLnPrior() would be used instead to set `last_ln_prior', but for moves that need
|	to revert this function is used because the previous log prior was saved and can thus easily be restored.
*/
void MCMCChainManager::setLastLnPrior(double ln_prior) 
	{
	last_ln_prior = ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Used by an MCMCUpdater-derived object to set the current log likelihood just before it returns from its update()
|	function.
*/
void MCMCChainManager::setLastLnLike(double ln_like) 
	{
	last_ln_like = ln_like;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current log joint prior. Each call to the calcJointLnPrior() member function refreshes the 
|	`last_ln_prior' data member.
*/
double MCMCChainManager::getLastLnPrior() const 
	{
	return last_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears all vectors (all_updaters, moves, model_params, edge_len_params and edge_len_hyperparams).
*/
void MCMCChainManager::clear()
	{
	all_updaters.clear();
	moves.clear();
	model_params.clear();
	edge_len_params.clear();
	edge_len_hyperparams.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `chain_mgr' data member of all parameters so that they can call MCMCChainManager::getLnPrior when any one 
|	of them needs to compute the joint prior over all parameters (or wishes to refresh the MCMCChainManager's 
|	`last_ln_like' data member. The MCMCUpdater::`chain_mgr' data member is a weak_ptr (if a shared_ptr were used, we 
|	would have a cycle: MCMCChainManager holds a shared_ptr to the MCMCUpdater object, and MCMCUpdater would hold a 
|	shared_ptr to MCMCChainManager). The weak_ptr that is passed to the setChainManager member function of each 
|	MCMCUpdater object is created from a shared_ptr to this, which we can get because MCMCChainManager was derived from
|	boost::enable_shared_from_this<MCMCChainManager>. Using the shared_from_this() function presumes that there exists 
|	a shared_ptr already managing this. Thus, a MCMCChainManager object should not be created within C++ without using 
|	shared_ptr to manage it.
*/
void MCMCChainManager::finalize()
	{
	// Determine how big the all_updaters vector needs to be
	unsigned sz = (unsigned)(moves.size() + edge_len_params.size() + model_params.size() + edge_len_hyperparams.size());
	if (sz == 0)
		{
		throw XLikelihood("must add at least one updater object to the chain manager before calling finalize()");
		}

	// Clear all_updaters vector
	all_updaters.clear();
	all_updaters.resize(sz);

	// Add the moves first
	moves_end = std::copy(moves.begin(), moves.end(), all_updaters.begin());
	moves_begin = all_updaters.begin();

	// Next add the normal edge length parameters
	edgelens_end = std::copy(edge_len_params.begin(), edge_len_params.end(), moves_end);
	edgelens_begin = moves_end;     

	// Next add the edge length hyperparameters (if there are any)
	hyperparams_end = std::copy(edge_len_hyperparams.begin(), edge_len_hyperparams.end(), edgelens_end);
	hyperparams_begin = edgelens_end;

	// Next add the remaining model parameters
	model_params_end = std::copy(model_params.begin(), model_params.end(), hyperparams_end);
	model_params_begin = hyperparams_end;

	params_begin = edgelens_begin;
	params_end = model_params_end;

	// Call each updater's setChainManager member function, supplying a shared pointer to this object
	// This allows each updater to ask the MCMCChainManager to calculate the joint prior over all
	// parameters when needed for computing the posterior density during slice sampling
	std::for_each(all_updaters.begin(), all_updaters.end(), 
		boost::lambda::bind(&MCMCUpdater::setChainManager, *boost::lambda::_1, ChainManagerWkPtr(shared_from_this())));

	// Call each parameter's setCurrValueFromModel member function to make sure that their 
    // curr_value data members are up to date
	std::for_each(params_begin, params_end, 
		boost::lambda::bind(&MCMCUpdater::setCurrValueFromModel, *boost::lambda::_1));

	// Call each parameter's recalcPrior member function to make sure their curr_ln_prior data members are
	// up to date
	std::for_each(params_begin, params_end, 
		boost::lambda::bind(&MCMCUpdater::recalcPrior, *boost::lambda::_1));

	//@POL Wondering why we need all these vectors? It is convenient to have one vector of updaters, but not all of 
	// these updaters are necessarily created in MCMCChainManager::addMCMCUpdaters. When finalize() is called, we 
	// can be sure that all the updaters that are going to be defined have been defined, and at that point we can 
	// quickly build the vector all_updaters. In the future, if it turns out we never need moves, model_params, or
	// edge_len_params after this point, we can just clear them at this point because we have iterators now set to
	// the beginning and end of each section of the all_updaters vector anyway.

	dirty = false;
	}

// Using USE_PARTIAL_SUM_FUNCTOR is slightly faster, but makes the assumption that the edge length prior
// distribution is Exponential, which is not good. This assumption is not necessary, but without it
// speed suffers and the alternative is preferable.
// 
// Note: now that separate priors have been installed for internal and external edges, the PartialPriorSum 
// approach is even less useful and has been abandoned. I've left the code here in case in case something
// like this needs to be done elsewhere in the future.
//#define USE_PARTIAL_SUM_FUNCTOR

//#if defined(USE_PARTIAL_SUM_FUNCTOR)
//struct PartialPriorSum : public std::unary_function<double, void>
//	{
//	double sum;
//	ExponentialDistribution & distr;
//
//	PartialPriorSum(ExponentialDistribution & p) : sum(0.0), distr(p) 
//		{}
//
//	void operator()(const double & nxt)
//		{
//	try 
//		{
//		sum += distr.GetLnPDF(nxt);
//		}
//	catch(XProbDist & x)
//		{
//		sum += distr.GetRelativeLnPDF(nxt);
//		}
//		}
//
//	double result()
//		{
//		return sum;
//		}
//	};
//#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the unnormalized prior for the supplied edge length using the prior for external edge lengths.
*/
double MCMCChainManager::calcExternalEdgeLenPriorUnnorm(
  double v) const
	{
    // The first edge length parameter always governs external edge lengths
    PHYCAS_ASSERT((edgelens_end - edgelens_begin) >= 1);
	MCMCUpdaterShPtr u = *edgelens_begin;
	ProbDistShPtr p = u->getPriorDist();
	double tmp = 0.0;
	try 
		{
		tmp = p->GetLnPDF(v);
		}
	catch(XProbDist &)
		{
		tmp =  p->GetRelativeLnPDF(v);
		}
    return tmp;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the unnormalized prior for the supplied edge length using the prior for internal edge lengths.
*/
double MCMCChainManager::calcInternalEdgeLenPriorUnnorm(double v) const
	{
    // The first edge length parameter governs internal edge lengths if there is only one edge length parameter.
    // If there are two edge length parameters, then the second one always governs internal edge lengths
	bool two_edgelen_priors = (edgelens_end - edgelens_begin == 2);
    MCMCUpdaterShPtr u = (two_edgelen_priors ? (*(edgelens_begin + 1)) : (*edgelens_begin));
	ProbDistShPtr p = u->getPriorDist();
	double tmp = 0.0;
	try 
		{
		tmp = p->GetLnPDF(v);
		}
	catch(XProbDist &)
		{
		tmp =  p->GetRelativeLnPDF(v);
		}
    return tmp;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Visits each parameter and calls its getLnPrior function to compute the joint log prior over all parameters. The
|	accumulate algorithm requires its fourth argument to be a binary functor, for which an unnamed boost lambda function
|	object is used. This function refreshes the value of the data member `last_ln_prior'.
*/
double MCMCChainManager::calcJointLnPrior()
	{
	if (all_updaters.empty())
		{
		throw XLikelihood("should not call MCMCChainManager::calcJointLnPrior before calling MCMCChainManager::finalize");
		}

	last_ln_prior = std::accumulate(params_begin, params_end, 0.0, 
		boost::lambda::_1 += boost::lambda::bind(&MCMCUpdater::getLnPrior, *boost::lambda::_2));

	return last_ln_prior;
	}

#if 0
/*----------------------------------------------------------------------------------------------------------------------
|	Uses the two edge length parameters in the `all_updaters' vector to compute the log-prior for each value in the 
|	supplied vector of edge lengths, returning the sum. Useful for efficiency purposes for moves that modify only a few 
|	edge lengths.
*/
double MCMCChainManager::partialEdgeLenPrior(const std::vector<double> & edge_len_vect) const
	{
    // To do next: replace all calls to this function (partialEdgeLenPrior) with calls to new functions
    // named calcInternalEdgeLenPriorUnnorm and calcExternalEdgeLenPriorUnnorm. There are only a couple of places
    // where the edge_len_vect vector has more than one element, so this function is not saving much
    // and we don't know which of the edge lengths supplied belong to internal or external edges.

	if (edgelens_end - edgelens_begin == 0)
		{
		throw XLikelihood("partialEdgeLenPrior requires a master edge length prior to have been added to the chain manager");
		}

	MCMCUpdaterShPtr u = *edgelens_begin;     
	ProbDistShPtr p = u->getPriorDist();

//#if defined(USE_PARTIAL_SUM_FUNCTOR)
//	ExponentialDistribution * d = dynamic_cast<ExponentialDistribution *>(p.get());
//	PHYCAS_ASSERT(d);
//	//@POL assuming edge length distribution is Exponential for now (allows GetLnPDF to be inlined)
//	return std::for_each(edge_len_vect.begin(), edge_len_vect.end(), PartialPriorSum(*d)).result();
//#else
	double partial_prior_sum = 0.0;
	for (std::vector<double>::const_iterator it = edge_len_vect.begin(); it != edge_len_vect.end(); ++it)
		{
		double tmp = 0.0;
		try 
			{
			tmp = p->GetLnPDF(*it);
			}
		catch(XProbDist &)
			{
			tmp =  p->GetRelativeLnPDF(*it);
			}
		partial_prior_sum += tmp;
		}
	return partial_prior_sum;
//#endif
	}
#endif

}	// namespace phycas
