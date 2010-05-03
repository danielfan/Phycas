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

#include <fstream>//temp 

#include <iterator>		// for std::distance 
#include <algorithm>	// for std::lower_bound  

//#include "phycas/force_include.h"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/mapping_move.hpp"

namespace phycas
{

void MCMCChainManager::debugUpdaterReport(std::string s)
	{
	std::cerr << "\n>>>>> MCMCChainManager updater report (" << s << ") <<<<<" << std::endl;
	std::cerr << "all_updaters:" << std::endl;
	for (MCMCUpdaterIter iter = all_updaters.begin(); iter != all_updaters.end(); ++iter)
		{
		std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
		}
	std::cerr << "moves:" << std::endl;
	for (MCMCUpdaterIter iter = moves.begin(); iter != moves.end(); ++iter)
		{
		std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
		}
	std::cerr << "model_params:" << std::endl;
	for (MCMCUpdaterIter iter = model_params.begin(); iter != model_params.end(); ++iter)
		{
		std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
		}
	std::cerr << "edge_len_params:" << std::endl;
	for (MCMCUpdaterIter iter = edge_len_params.begin(); iter != edge_len_params.end(); ++iter)
		{
		std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
		}
	std::cerr << "edge_len_hyperparams:" << std::endl;
	for (MCMCUpdaterIter iter = edge_len_hyperparams.begin(); iter != edge_len_hyperparams.end(); ++iter)
		{
		std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
		}
	std::cerr << "\n" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor for 
*/
MCMCChainManager::~MCMCChainManager() 
	{
	//std::cerr << "\n>>>>> MCMCChainManager dying..." << std::endl;
	releaseUpdaters();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears vectors of shared pointers to MCMCUpdater objects. Called from destructor (or can be called at will). 
*/
void MCMCChainManager::releaseUpdaters()
    {
	//std::cerr << "\n>>>>> MCMCChainManager::releaseUpdaters()..." << std::endl;    
	//std::cerr << "\nBefore all_updaters cleared..." << std::endl;
	for (MCMCUpdaterIter iter = all_updaters.begin(); iter != all_updaters.end(); ++iter)
		{
		//std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
		(*iter)->releaseSharedPointers();
		}
	all_updaters.clear();

	//std::cerr << "\nBefore moves cleared..." << std::endl;
	//for (MCMCUpdaterIter iter = moves.begin(); iter != moves.end(); ++iter)
	//	{
	//	std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
	//	}
	moves.clear();

	//std::cerr << "\nBefore model_params cleared..." << std::endl;
	//for (MCMCUpdaterIter iter = model_params.begin(); iter != model_params.end(); ++iter)
	//	{
	//	std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
	//	}
	model_params.clear();

	//std::cerr << "\nBefore edge_len_params cleared..." << std::endl;
	//for (MCMCUpdaterIter iter = edge_len_params.begin(); iter != edge_len_params.end(); ++iter)
	//	{
	//	std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
	//	}
	edge_len_params.clear();

	//std::cerr << "\nBefore edge_len_hyperparams cleared..." << std::endl;
	//for (MCMCUpdaterIter iter = edge_len_hyperparams.begin(); iter != edge_len_hyperparams.end(); ++iter)
	//	{
	//	std::cerr << "  " << (*iter)->getName() << " use count = " << iter->use_count() << std::endl;
	//	}
	edge_len_hyperparams.clear();
    }


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

    //@POL should only go to this trouble if univents are being used
    // Try to find a MappingMove updater to use for this. If one is present, its update method must
    // be called before the likelihood can be calculated in order to generate the latent variables used in
    // uniformized mapping MCMC
    bool found = false; 
	for (MCMCUpdaterIter it = moves_begin; it != moves_end; ++it)
        {
        MappingMove * m = dynamic_cast<MappingMove *>(it->get());
        if (m)
            {
            m->update();
            last_ln_like = m->recalcLike();
            found = true;
            break;
            }
        }
    if (!found)
        {
        MCMCUpdaterShPtr u = *params_begin;
	    last_ln_like = u->recalcLike();
        }
#if 0 //@POL this is what should be executed if not doing unimap 
	MCMCUpdaterShPtr u = *params_begin;
	last_ln_like = u->recalcLike();
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates and returns the log of the joint working prior by calling recalcWorkingPrior() for all parameters.
*/
double MCMCChainManager::recalcLnWorkingPrior() const
	{
	// Might need to add more checks here
	if (dirty)
		{
		throw XLikelihood("cannot call recalcLnWorkingPrior() for chain manager before calling finalize()");
		}
	
	//std::cerr << "\nIn MCMCChainManager::recalcLnWorkingPrior()..." << std::endl;//temp
	
	// Visit all updaters and let those who are prior stewards update their component of the joint prior.
	double ln_working_prior = 0.0;
	for (MCMCUpdaterConstIter it = all_updaters.begin(); it != all_updaters.end(); ++it)
		{
		const boost::shared_ptr<MCMCUpdater> s = *it;
		if (s->isPriorSteward())
			{
			double this_ln_prior = s->recalcWorkingPrior();
			ln_working_prior += this_ln_prior;
			//std::cerr << boost::str(boost::format("@@@@@@@@@@ this working prior = %g, cumulative = %g, name = '%s'") % this_ln_prior % ln_working_prior % s->getName()) << std::endl;
			}
		}
	return ln_working_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Refreshes `last_ln_prior' by calling recalcPrior() for all parameters.
*/
void MCMCChainManager::refreshLastLnPrior()
	{
	last_ln_prior = 0.0;
	if (samc_likelihood_only)
		{
		return;
		}
		
	// Might need to add more checks here
	if (dirty)
		{
		throw XLikelihood("cannot call refreshLastLnPrior() for chain manager before calling finalize()");
		}
		
	// Visit all updaters and let those who are prior stewards update their component of the joint prior.
	for (MCMCUpdaterConstIter it = all_updaters.begin(); it != all_updaters.end(); ++it)
		{
		const boost::shared_ptr<MCMCUpdater> s = *it;
		if (s->isPriorSteward() && !s->isFixed())
			{
			double this_ln_prior = s->recalcPrior();
			last_ln_prior += this_ln_prior;
			//std::cerr << boost::str(boost::format("@@@@@@@@@@ Current prior = %g, Cumulative = %g, name = '%s'") % this_ln_prior % last_ln_prior % s->getName()) << std::endl;
			}
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
  unsigned weight,					/**< is the weight to be used for all parameters added by this function */
  int subset_pos)					/**< if 0 (first subset) or -1 (only subset), edge length parameters and hyperparams will be added */
	{
	//std::cerr << "\n***** model.getExternalEdgeLenPrior() upon entering MCMCChainManager::addMCMCUpdaters: use count = " << m->getExternalEdgeLenPrior().use_count() << std::endl;
	
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
	if (subset_pos <= 0)
		{
		if (!m->getInternalEdgeLenPrior())
			{
			throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no internal edge length prior defined");
			}
		if (!m->getExternalEdgeLenPrior())
			{
			throw XLikelihood("Error in MCMCChainManager::addMCMCUpdaters: no external edge length prior defined");
			}
		}

	MCMCUpdaterIter iter;
	MCMCUpdaterVect edgelens;
	MCMCUpdaterVect edgelen_hyperparams;
	MCMCUpdaterVect parameters;

	// Ask the model to create the edge length parameters, edge length hyperparameters, and 
	// its own model-specific parameters.
	m->createParameters(t, edgelens, edgelen_hyperparams, parameters, subset_pos);

	if (subset_pos <= 0)
		{
		// Add the edge length parameters. If length 2, then first will be the master parameter for external edges
		// and the second will be the master parameter for the external edges. If length 1, then it will be the
		// edge length master parameter that governs all edges in the tree.
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

		// Add the edge length hyperparameters (if any were created). If length 2, then first will be the hyperparameter 
		// governing the mean of the external edge length prior and the second will be the hyperparameter governing
		// the mean of the internal edge length prior. If length 1, then it will be the hyperparameter governing the
		// mean of the prior for all edge lengths.
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
		}

	//std::cerr << boost::str(boost::format("~~~ length of parameters vector is %d") % parameters.size()) << std::endl;

	for (iter = parameters.begin(); iter != parameters.end(); ++iter)
		{
		(*iter)->setWeight(weight);
		(*iter)->setMaxUnits(max_units);
		(*iter)->setModel(m);
		(*iter)->setTreeLikelihood(like);
		(*iter)->setLot(r);
		//std::cerr << boost::str(boost::format("~~~ calling MCMCChainManager::addModelParam for parameter %s") % (*iter)->getName()) << std::endl;
		addModelParam(*iter);
		}
	}

////////////////////////////////////////// inserted contents of mcmc_chain_manager.inl here /////////////////////////////

/*----------------------------------------------------------------------------------------------------------------------
|	The MCMCChainManager constructor does nothing.
*/
MCMCChainManager::MCMCChainManager() 
  : last_ln_like(0.0)
  , last_ln_prior(0.0)
  , doing_samc(false)
  , samc_likelihood_only(false)
  , samc_best(-DBL_MAX)
  , samc_rfbest(UINT_MAX)
  , samc_gain(1.0)
  , samc_t0(50.0)
  , samc_eta(0.6)
  , samc_proposal_matrix(0)
  , dirty(true)
	{
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `doing_samc' which is set to true if initSAMC is called with a non-empty vector
|	of log energy levels and is set to false if initSAMC is never called or is called with an empty vector.
*/
bool MCMCChainManager::doingSAMC() const
	{
	return doing_samc;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Cleans up SAMC-related data structures.
*/
void MCMCChainManager::finalizeSAMC()
	{
	std::cerr << ">>>>>>>>> In MCMCChainManager::finalizeSAMC <<<<<<<<<<" << std::endl;
	if (samc_proposal_matrix)
		{
		unsigned i, j;
		unsigned sz = (unsigned)samc_count.size();
		
		// save the proposal matrix to a file
		std::cerr << ">>>>>>>>> Saving the SAMC proposal matrix to file proposal_matrix.txt <<<<<<<<<<" << std::endl;
		std::ofstream doof("proposal_matrix.txt");
		for (i = 0; i < sz; ++i)
			{
			for (j = 0; j < sz; ++j)
				{
				doof << boost::str(boost::format("%12d\t") % samc_proposal_matrix[i][j]);
				}			
			doof << std::endl;
			}
		doof.close();
		
		// delete storage for the proposal matrix
		for (i = 0; i < sz; ++i)
			{
			delete [] samc_proposal_matrix[i];
			}
		delete [] samc_proposal_matrix;
		samc_proposal_matrix = 0;
		}
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `energy_levels' data member to a copy of the supplied values. If `elevels' is empty, clears `energy_levels'.
*/
void MCMCChainManager::initSAMC(
  const double_vect_t & elevels,
  double t0,
  double eta)
	{
	if (elevels.empty())
		{
		doing_samc = false;
		samc_loglevels.clear();
		samc_pi.clear();
		samc_theta.clear();
		samc_count.clear();
		samc_t0 = 50.0;
		samc_eta = 0.6;
		}
	else
		{
		doing_samc = true;
		samc_t0 = t0;
		samc_eta = eta;
		unsigned m = (unsigned)(elevels.size() + 1);
		double freq = 1.0/(double)m;
		
		samc_loglevels.resize(m - 1);
		std::copy(elevels.begin(), elevels.end(), samc_loglevels.begin());
		
		std::cerr << "\n>>>>> elevels: ";
		std::copy(elevels.begin(), elevels.end(), std::ostream_iterator<double>(std::cerr, " "));
		std::cerr << "\n>>>>> loglevels: ";
		std::copy(samc_loglevels.begin(), samc_loglevels.end(), std::ostream_iterator<double>(std::cerr, " "));
		std::cerr << "\n" << std::endl;
	
		samc_pi.resize(m);
#if 0
		samc_pi.assign(m, freq);
#else
		double factor = 2.0;
		double largest_factor = exp(log(DBL_MAX)/(double)(m + 1));
		if (largest_factor > factor)
			factor = largest_factor;
		std::cerr << "@#@#@#@#@#@#@#@#@#@#@# m              = " << m << " @#@#@#@#@#@#@#@#@#@#@#" << std::endl;
		std::cerr << "@#@#@#@#@#@#@#@#@#@#@# factor         = " << factor << " @#@#@#@#@#@#@#@#@#@#@#" << std::endl;
		std::cerr << "@#@#@#@#@#@#@#@#@#@#@# largest_factor = " << largest_factor << " @#@#@#@#@#@#@#@#@#@#@#" << std::endl;
		double norm_const = (pow(factor, (double)m) - 1.0)/(factor - 1.0);
		std::cerr << "@#@#@#@#@#@#@#@#@#@#@# norm_const = " << norm_const << " @#@#@#@#@#@#@#@#@#@#@#" << std::endl;
		double curr = 1.0;
		double check_sum = 0.0;
		for (unsigned i = 0; i < m; ++i)
			{
			samc_pi[i] = curr/norm_const;
			std::cerr << "@#@#@#@#@#@#@#@#@#@#@# " << boost::str(boost::format("samc_pi[%d] = %.9f") % i % samc_pi[i]) << " @#@#@#@#@#@#@#@#@#@#@#" << std::endl;
			curr *= factor;
			check_sum += samc_pi[i];
			}
		PHYCAS_ASSERT(fabs(check_sum - 1.0) < 0.001);
		std::cerr << "@#@#@#@#@#@#@#@#@#@#@# check_sum      = " << check_sum << " @#@#@#@#@#@#@#@#@#@#@#" << std::endl;
#endif
		
		samc_theta.resize(m);
		samc_theta.assign(m, 0.0);
		
		samc_count.resize(m);
		samc_count.assign(m, 0);
		}
	}
		
/*----------------------------------------------------------------------------------------------------------------------
|	Adds 1 to the count of proposals from level `from_level' to level `to_level'.
*/
void MCMCChainManager::debugAddSAMCProposalPair(
  unsigned from_level,	/**< is the current energy level */
  unsigned to_level)	/**< is the proposed level */
	{
	if (!samc_proposal_matrix)
		{
		unsigned sz = (unsigned)samc_count.size();
		samc_proposal_matrix = new unsigned *[sz];
		for (unsigned i = 0; i < sz; ++i)
			{
			samc_proposal_matrix[i] = new unsigned[sz];
			for (unsigned j = 0; j < sz; ++j)
				samc_proposal_matrix[i][j] = 0;
			}
		}
	++samc_proposal_matrix[from_level][to_level];
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the index of the energy level of the most recent sampled point based on the supplied `logf'.
*/
unsigned MCMCChainManager::getSAMCEnergyLevel(
  double logf)	const 	/**< is the log posterior of the current point */
	{
	PHYCAS_ASSERT(doing_samc);
	
	// Suppose logf = -0.51 and supposing also that these are the n = 5 levels:
	//
	// level 4
	// --------  0.00
	// level 3
	// -------- -0.29
	// level 2        <-- -0.51
	// -------- -0.69
	// level 1
	// -------- -1.39
	// level 0
	//
	// then x_index = 2.
	double_vect_t::const_iterator it = std::lower_bound(samc_loglevels.begin(), samc_loglevels.end(), logf);
	unsigned x_index = (unsigned)samc_loglevels.size();
	if (it != samc_loglevels.end())
		x_index = (unsigned)std::distance(samc_loglevels.begin(), it);
	//else 
	//	{
	//	std::cerr << "\n>>>>> loglevels: ";
	//	std::copy(samc_loglevels.begin(), samc_loglevels.end(), std::ostream_iterator<double>(std::cerr, " "));
	//	std::cerr << "\n>>>>> logf    = " << logf;
	//	std::cerr << "\n>>>>> x_index = " << x_index;
	//	std::cerr << "\n" << std::endl;
	//	}
	
	return x_index;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the SAMC weight of the most recent sampled point based on the supplied energy level index `i'.
*/
double MCMCChainManager::getSAMCWeight(
  unsigned i) const 	/**< is the index of the energy level of the current point */
	{
	PHYCAS_ASSERT(doing_samc);
	PHYCAS_ASSERT(i < samc_theta.size());
	return samc_theta[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `samc_gain'.
*/
double MCMCChainManager::getSAMCGain() const
	{
	PHYCAS_ASSERT(doing_samc);
	return samc_gain;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `samc_best'.
*/
double MCMCChainManager::getSAMCBest() const
	{
	PHYCAS_ASSERT(doing_samc);
	return samc_best;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the vector of SAMC weights (data member `samc_theta') associated with each SAMC energy level.
*/
double_vect_t MCMCChainManager::getSAMCWeights() const
	{
	PHYCAS_ASSERT(doing_samc);
	return samc_theta;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the vector of SAMC log energy level boundaries (data member `samc_loglevels').
*/
double_vect_t MCMCChainManager::getSAMCLogLevels() const
	{
	PHYCAS_ASSERT(doing_samc);
	return samc_loglevels;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the vector of SAMC energy level proportions (data member `samc_pi') associated with each SAMC energy level.
*/
double_vect_t MCMCChainManager::getSAMCPi() const
	{
	PHYCAS_ASSERT(doing_samc);
	return samc_pi;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the vector of SAMC energy level proportions (data member `samc_pi') associated with each SAMC energy level.
*/
unsigned MCMCChainManager::calcRFDistance(TreeShPtr ref_tree) const
	{
	PHYCAS_ASSERT(doing_samc);
	MCMCUpdaterShPtr u = *params_begin;
	return u->getTree()->robinsonFoulds(ref_tree);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `samc_rfbest', which is the last-calculated Robinson-Foulds distance
|	between the SAMC reference tree and the best tree seen thus far.
*/
unsigned MCMCChainManager::getSAMCRobinsonFouldsBest() const
	{
	PHYCAS_ASSERT(doing_samc);
	return samc_rfbest;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `samc_likelihood_only', which determines whether the normalized
|	likelihood or the posterior will be sampled during SAMC analyses.
*/
bool MCMCChainManager::getSAMCLikelihoodOnly() const
	{
	return samc_likelihood_only;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the root mean square error of the energy level proportions. This serves as a measure of how well the SAMC
|	analysis is converging.
*/
double MCMCChainManager::getSAMCPiRMSE() const
	{
	PHYCAS_ASSERT(doing_samc);
	
	double rmse = 0.0;
	unsigned nlevels = (unsigned)samc_count.size();
	PHYCAS_ASSERT(nlevels > 0);
	unsigned total_count = std::accumulate(samc_count.begin(), samc_count.end(), 0);
	for (unsigned i = 0; i < nlevels; ++i)
		{
		double term = ((double)samc_count[i]/(double)total_count) - samc_pi[i];
		rmse += pow(term, 2.0);
		//std::cerr << boost::str(boost::format("term = %g, count = %d, pi = %g, rmse = %g") % term % samc_count[i] % samc_pi[i] % rmse) << std::endl;
		}
	//std::cerr << std::endl;
	rmse /= (double)nlevels;
	return sqrt(rmse);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the counts of visits to each SAMC energy level.
*/
uint_vect_t MCMCChainManager::getSAMCCounts() const
	{
	PHYCAS_ASSERT(doing_samc);
	return samc_count;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of `samc_likelihood_only', which indicates, if true, that SAMC analyses should explore the likelihood
|	surface, not the posterior. If false, SAMC analyses will explore the usual posterior distribution.
*/
void MCMCChainManager::setSAMCLikelihoodOnly(
  bool likelihood_only_status) 	/**< is the desired setting for `samc_likelihood_only'  */
	{
	//temp! PHYCAS_ASSERT(doing_samc);
	samc_likelihood_only = likelihood_only_status;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of `samc_ref_tree', which should be a shared pointer to the best tree known.
*/
void MCMCChainManager::setSAMCRefTree(
  TreeShPtr ref_tree) 	/**< is the SAMC reference tree */
	{
	PHYCAS_ASSERT(doing_samc);
	samc_ref_tree = ref_tree;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of `samc_best', which should be the log posterior of the best point sampled by SAMC thus far.
*/
void MCMCChainManager::setSAMCBest(
  double new_best) 	/**< is the new best log posterior value */
	{
	PHYCAS_ASSERT(doing_samc);
	PHYCAS_ASSERT(new_best > samc_best);
	samc_best = new_best;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Raises (if `incr' is positive) or lowers (if `incr' is negative) each level boundary by an amount `incr'.
*/
void MCMCChainManager::adjustSAMCLevels(
  double incr,			/**< is the amount by which to adjust all level boundaries (if fix_levels_now is False) or the point at which to fix the highest level (if fix_levels_now is True) */
  bool fix_levels_now)	/**< determines whether we continue to raise the sea level or fix the energy level boundaries */
	{
	if (fix_levels_now)
		{
		double lowest  = *(samc_loglevels.begin());
		double highest = lowest + incr;
		
		// reset level boundaries
		samc_loglevels.clear();
		for (double logf = lowest; logf < highest; logf += 10.0)
			{
			samc_loglevels.push_back(logf);
			}
		unsigned nboundaries = (unsigned)samc_loglevels.size();
		unsigned nlevels = nboundaries + 1;
		PHYCAS_ASSERT(nlevels > 1);
		
		// reset counts
		samc_count.resize(nlevels);
		samc_count.assign(nlevels, 0);
		
		// reset weights
		samc_theta.resize(nlevels);
		samc_theta.assign(nlevels, 0.0);
		
		// reset stationary frequencies
		samc_pi.resize(nlevels);
		samc_pi.assign(nlevels, 1.0/(double)nlevels);
		}
	else 
		{
		std::transform(samc_loglevels.begin(), samc_loglevels.end(), samc_loglevels.begin(), boost::lambda::_1 + incr);
		unsigned sz = (unsigned)samc_count.size();
		samc_count.assign(sz, 0);
		}
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Updates the SAMC weights based on the supplied energy level `i' of the most recent sampled point.
*/
void MCMCChainManager::updateSAMCWeights(
  unsigned i)	/**< is the index of the level to which the current point belongs, with i == 0 representing the highest energy (lowest posterior density) */
	{
	PHYCAS_ASSERT(doing_samc);
	PHYCAS_ASSERT(i < samc_count.size());
	
	// update counts
	samc_count[i] += 1;
	
	// update weights        
	unsigned m = (unsigned)samc_theta.size();
	for (unsigned j = 1; j < m; ++j)
		{
		if (samc_count[j] > 0)
			{
			// The following leaves samc_theta[0] equal to its initial value, 0.0
			double i_equals_0 = (i == 0 ? 1.0 : 0.0);
			double i_equals_j = (i == j ? 1.0 : 0.0);
			double first_term  = i_equals_j - samc_pi[j];
			double second_term = i_equals_0 - samc_pi[0];
			samc_theta[j] += samc_gain*(first_term - second_term);
			}
		}
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Updates the SAMC gain based on the supplied time `t'.
*/
void MCMCChainManager::updateSAMCGain(
  unsigned t)	/**< is the current number of cycles */
	{
	PHYCAS_ASSERT(doing_samc);
	double denom = (t > samc_t0 ? t : samc_t0);
	double log_gain = samc_eta*(std::log(samc_t0) - std::log(denom));
	samc_gain = std::exp(log_gain);	
	}
#endif
	
#define DEBUG_UPDATERS  0
#define DEBUG_SAMC      0

#if DEBUG_SAMC
unsigned samc_trace_count = 0;
unsigned samc_trace_mod = 100;
#endif	
/*----------------------------------------------------------------------------------------------------------------------
|	For all updaters stored in `all_updaters', obtain the weight w and call the update fuction of the updater w times.
*/
void MCMCChainManager::updateAllUpdaters()
	{
#if DEBUG_SAMC
	std::ofstream tmpf("samc_trace.txt", std::ios::out | std::ios::app);
#endif	
	for (MCMCUpdaterVect::iterator it = all_updaters.begin(); it != all_updaters.end(); ++it)
		{
#if DEBUG_UPDATERS
		std::string nm = (*it)->getName();
		unsigned w = (*it)->getWeight();
		for (unsigned i = 0; i < w; ++i)
			{
			std::cerr << "########## updating " << nm << "..." << std::endl;
			(*it)->setSaveDebugInfo(true);
			(*it)->update();
			std::cerr << boost::str(boost::format("%s | %s") % nm % (*it)->getDebugInfo()) << std::endl;
			}
#else
		unsigned w = (*it)->getWeight();
		for (unsigned i = 0; i < w; ++i)
			{
			(*it)->update();
			
#if DEBUG_SAMC
			double tmp_last_ln_prior = getLastLnPrior();//temp!
			calcJointLnPrior();//temp!
			if (fabs(tmp_last_ln_prior - getLastLnPrior()) > 0.01)
				{
				std::cerr << "@@@@@@@@@@@@@@@@@@@@@@@ BAD BAD BAD BAD BAD BAD @@@@@@@@@@@@@@@@@@@@@@" << std::endl;
				}
#endif	
			double logf = getLastLnLike() + getLastLnPrior();
#if DEBUG_SAMC
			if (samc_trace_count % samc_trace_mod == 0)
				{
				tmpf << boost::str(boost::format("%.6f\n") % logf);
				samc_trace_count = 0;
				}
			else {
				++samc_trace_count;
			}

#endif	
			if (doing_samc)
				{
				if (logf > samc_best)
					{
					samc_best = logf;
					TreeShPtr curr_tree = (*it)->getTree();
					std::ofstream doof("samc_best.tre");
					doof << "#nexus\n\nbegin trees;\n  tree samc_best = [&U] " << curr_tree->MakeNewick() << ";\nend;" << std::endl;
					doof.close();
					}
				if (samc_ref_tree)
					{
					TreeShPtr curr_tree = (*it)->getTree();
					unsigned tmp = curr_tree->robinsonFoulds(samc_ref_tree);
					if (tmp == 0 && samc_rfbest > 0)
						{
						std::ofstream doof("samc_first_rf0.tre");
						doof << "#nexus\n\nbegin trees;\n  tree samc_first_rf0 = [&U] " << curr_tree->MakeNewick() << ";\nend;" << std::endl;
						doof.close();
						}
					if (tmp < samc_rfbest)
						{
						samc_rfbest = tmp;
						std::ofstream doof("samc_best_rf.tre");
						doof << "#nexus\n\nbegin trees;\n  tree samc_best_rf = [&U] " << curr_tree->MakeNewick() << ";\nend;" << std::endl;
						doof.close();
						}
					}
				}
			} 
#endif
		}
#if DEBUG_SAMC
	tmpf.close();
#endif	
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
void MCMCChainManager::addEdgeLenHyperparam(MCMCUpdaterShPtr p) 
	{
	edge_len_hyperparams.push_back(p);
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the specified value in the `edge_len_hyperparams' vector. If there are two edge length hyperparameters, the 
|   first applies to external edge lengths and the second to internal edge lengths.
*/
void MCMCChainManager::setEdgeLenHyperparam(
  unsigned i,   /**< is the index of the value to set */
  double v)     /**< is the value to set */
	{
    PHYCAS_ASSERT(i < edge_len_hyperparams.size());
	edge_len_hyperparams[i]->sendCurrValueToModel(v);
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
	//std::for_each(all_updaters.begin(), all_updaters.end(), 
	//	boost::lambda::bind(&MCMCUpdater::setChainManager, *boost::lambda::_1, ChainManagerWkPtr(shared_from_this())));
	ChainManagerWkPtr wptr(shared_from_this());
	for (MCMCUpdaterIter uit = all_updaters.begin(); uit != all_updaters.end(); ++uit)
		{
		//std::string nm = (*uit)->getName();
		//unsigned uc1 = (unsigned)(*uit).use_count();
		(*uit)->setChainManager(wptr);
		//unsigned uc2 = (unsigned)(*uit).use_count();
		//std::cerr << ">>>> updater " << nm << " (" << uc1 << " --> " << uc2 << ")" << std::endl;
		}
	
	for (MCMCUpdaterIter piter = params_begin; piter != params_end; ++piter)
		{
		// Call each parameter's recalcPrior member function to make sure their curr_ln_prior data members are
		// up to date
		(*piter)->recalcPrior();
		}

	//@POL Wondering why we need all these vectors? It is convenient to have one vector of updaters, but not all of 
	// these updaters are necessarily created in MCMCChainManager::addMCMCUpdaters. When finalize() is called, we 
	// can be sure that all the updaters that are going to be defined have been defined, and at that point we can 
	// quickly build the vector all_updaters. In the future, if it turns out we never need moves, model_params, or
	// edge_len_params after this point, we can just clear them at this point because we have iterators now set to
	// the beginning and end of each section of the all_updaters vector anyway.

	dirty = false;
	
#if 0
	std::cerr << std::endl;
	std::cerr << "##############################################" << std::endl;
	std::cerr << "Here are the updaters that can compute priors:" << std::endl;
	for (MCMCUpdaterIter it = all_updaters.begin(); it != all_updaters.end(); ++it)
		{
		MCMCUpdaterShPtr s = *it;
		if (s->computesUnivariatePrior())
			std::cerr << boost::str(boost::format("  '%s' (univariate prior)") % s->getName()) << std::endl;
		else if (s->computesMultivariatePrior())
			std::cerr << boost::str(boost::format("  '%s' (multivariate prior)") % s->getName()) << std::endl;
		else
			std::cerr << boost::str(boost::format("  '%s' (does not compute prior)") % s->getName()) << std::endl;
		}
		
	std::cerr << std::endl;
	std::cerr << "##############################################" << std::endl;
	std::cerr << "Here is the list of just the parameters:" << std::endl;
	for (MCMCUpdaterIter p = params_begin; p != params_end; ++p)
		{
		MCMCUpdaterShPtr s = *p;
		if (s->computesUnivariatePrior())
			std::cerr << boost::str(boost::format("  '%s' (univariate prior)") % s->getName()) << std::endl;
		else if (s->computesMultivariatePrior())
			std::cerr << boost::str(boost::format("  '%s' (multivariate prior)") % s->getName()) << std::endl;
		else
			std::cerr << boost::str(boost::format("  '%s' (does not compute prior)") % s->getName()) << std::endl;
		}
#endif
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
	if (samc_likelihood_only)
		return 0.0;
		
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
		PHYCAS_ASSERT(0);
		}
    return tmp;
    }

#if USING_EDGE_SPECIFIC_WORKING_PRIORS
/*----------------------------------------------------------------------------------------------------------------------
|	Computes the log of the working prior for the supplied external edge length.
*/
double MCMCChainManager::calcExternalEdgeLenWorkingPrior(
  const TreeNode & nd,	/**< is the node for which the working prior log density should be computed */
  double v) const		/**< is the edge length for which the working prior log density should be computed */
	{
	// If there are two edge_length_parameters, first is for external and second is for internal edges
	// (see MCMCChainManager::addMCMCUpdaters). If only one, it handles all edge lengths, so in this
	// case (external edge length) we know we should go with the first one.
    const MCMCUpdaterVect & edge_length_params = getEdgeLenParams();
	EdgeLenMasterParam * p = dynamic_cast<EdgeLenMasterParam *>(edge_length_params[0].get());
	PHYCAS_ASSERT(p);
	return p->lnWorkingPriorOneEdge(nd, v);
	}
#else
/*----------------------------------------------------------------------------------------------------------------------
|	Computes the log of the working prior for the supplied external edge length.
*/
double MCMCChainManager::calcExternalEdgeLenWorkingPrior(
  double v) const	/**< is the edge length for which the working prior log density should be computed */
	{
	// If there are two edge_length_parameters, first is for external and second is for internal edges
	// (see MCMCChainManager::addMCMCUpdaters). If only one, it handles all edge lengths, so in this
	// case (external edge length) we know we should go with the first one.
    const MCMCUpdaterVect & edge_length_params = getEdgeLenParams();
	return edge_length_params[0]->calcLnWorkingPriorPDF(v);
	}
#endif
	
#if USING_EDGE_SPECIFIC_WORKING_PRIORS
/*----------------------------------------------------------------------------------------------------------------------
|	Computes the log of the working prior for the supplied internal edge length.
*/
double MCMCChainManager::calcInternalEdgeLenWorkingPrior(
  const TreeNode & nd,		/**< is the node for which the working prior log density should be computed */
  double v) const			/**< is the edge length for which the working prior log density should be computed */
	{	
	// If there are two edge_length_parameters, first is for external and second is for internal edges
	// (see MCMCChainManager::addMCMCUpdaters). If only one, it handles all edge lengths.
	const MCMCUpdaterVect & edge_length_params = getEdgeLenParams();
	if (edge_length_params.size() == 2)
		{
		EdgeLenMasterParam * p = dynamic_cast<EdgeLenMasterParam *>(edge_length_params[1].get());
		PHYCAS_ASSERT(p);
		return p->lnWorkingPriorOneEdge(nd, v);
		//return edge_length_params[1]->lnWorkingPriorOneEdge(nd, false);
		}
	else 
		{
		EdgeLenMasterParam * p = dynamic_cast<EdgeLenMasterParam *>(edge_length_params[0].get());
		PHYCAS_ASSERT(p);
		return p->lnWorkingPriorOneEdge(nd, v);
		//return edge_length_params[0]->lnWorkingPriorOneEdge(nd, false);
		}
	}
#else
/*----------------------------------------------------------------------------------------------------------------------
|	Computes the log of the working prior for the supplied internal edge length.
*/
double MCMCChainManager::calcInternalEdgeLenWorkingPrior(
  double v) const	/**< is the edge length for which the working prior log density should be computed */
	{	
	// If there are two edge_length_parameters, first is for external and second is for internal edges
	// (see MCMCChainManager::addMCMCUpdaters). If only one, it handles all edge lengths.
	const MCMCUpdaterVect & edge_length_params = getEdgeLenParams();
	if (edge_length_params.size() == 2)
		return edge_length_params[1]->calcLnWorkingPriorPDF(v);
	else 
		return edge_length_params[0]->calcLnWorkingPriorPDF(v);
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the unnormalized prior for the supplied edge length using the prior for internal edge lengths.
*/
double MCMCChainManager::calcInternalEdgeLenPriorUnnorm(double v) const
	{
	if (samc_likelihood_only)
		return 0.0;
		
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
		PHYCAS_ASSERT(0);
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
	//	if (all_updaters.empty())
	//		{
	//		throw XLikelihood("should not call MCMCChainManager::calcJointLnPrior before calling MCMCChainManager::finalize");
	//		}
	//
	//	last_ln_prior = std::accumulate(params_begin, params_end, 0.0, 
	//		boost::lambda::_1 += boost::lambda::bind(&MCMCUpdater::getLnPrior, *boost::lambda::_2));
	refreshLastLnPrior();
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
		    PHYCAS_ASSERT(0);
			}
		partial_prior_sum += tmp;
		}
	return partial_prior_sum;
//#endif
	}
#endif

}	// namespace phycas
