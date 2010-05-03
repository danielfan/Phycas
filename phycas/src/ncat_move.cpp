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

#if defined(FLEXCAT_MODEL)

//#include "phycas/force_include.h"
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/ncat_move.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/tree_manip.hpp"

//#define SHOW_DEBUGGING_OUTPUT

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor calls NCatMove::reset to initialize all variables.
*/
NCatMove::NCatMove() : 
	MCMCUpdater(), 
	ncat_max(1),
	lambda(1.0),
	s(1),
	L(1.0),
	phi(0.5), 
	total_updates(0), 
  	ncat_distr(1, 0)
	{
	reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Re-initializes all data members that are recomputed each time ProposeNewState is called.
*/
void NCatMove::reset()
	{
	addcat_move_proposed	= false;
	u1						= 0.0;
	tmp_rate				= 0.0;
	tmp_prob				= 0.0;
	ln_jacobian				= 0.0;
	ln_hastings				= 0.0;
	ncat_before				= 0;
	ncat_after				= 0;
	rate_left				= 0.0;
	rate_right				= 0.0;
	//@POL tmp_rate_iter and tmp_prob_iter are not being reset - does that matter?
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Proposes either an addcat move with probability phi or a delcat move with probability 1 - phi. If there is only 1
|	rate category, proposes an addcat move with probability 1.
*/
void NCatMove::proposeNewState()
	{
#if DISABLED_UNTIL_FLEXCAT_WORKING_WITH_PARTITIONING
	ncat_before = model->getNGammaRates();
	addcat_move_proposed = (ncat_before == 1 ? true : (bool)(rng->Uniform(FILE_AND_LINE) < phi));

	// Make sure model is using the same upper bound on unnormalized rates that we are
	model->setFlexRateUpperBound(L);

	if (addcat_move_proposed)
		{
		ncat_after = ncat_before + 1;

		// set u1, tmp_rate and tmp_prob data members
		proposeAddCatMove();

		// Compute the log of the Jacobian
		ln_jacobian = std::log(L);

		// Compute the log of the Hastings ratio 
		ln_hastings = std::log(1.0 - phi) - std::log((double)ncat_after);
		if (ncat_before > 1)
			ln_hastings -= std::log(phi);
		}
	else
		{
		ncat_after = ncat_before - 1;
		PHYCAS_ASSERT(ncat_after > 0);

		proposeDelCatMove();

		// Compute the log of the Jacobian
		ln_jacobian = -std::log(L);

		// Compute the log of the Hastings ratio 
		ln_hastings = std::log((double)ncat_before) - std::log(1.0 - phi);
		if (ncat_after > 1)
			ln_hastings += std::log(phi);
		}

	// If number of rates has increased beyond the current capacity of CLAs, be sure to 
	// store all CLAs checked out to the tree now so that they will be deleted when  
	// likelihood->cla_pool.clearStack() is called (in the call to recalcRelativeRates below)
	unsigned new_nr = model->getNRatesTotal();
	const CondLikelihoodStorageShPtr clapool = likelihood->getCLAStorage();
	unsigned old_nr = clapool->getNumRates();

	if (new_nr > old_nr)
		{
		likelihood->storeAllCLAs(tree);
		}

	// Renormalize the rates and probs so that we are ready for next likelihood calculation
	likelihood->recalcRelativeRates();

	// If the number of rate categories has increased, check to see if it is now greater than
	// the maximum seen thus far. If so, need to reallocate transition probability matrices
	// in the tree to accommodate the larger number of rate categories
	if (ncat_after > ncat_max)
		{
		ncat_max = ncat_after;
		likelihood->prepareForLikelihood(tree);
		}
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void NCatMove::proposeAddCatMove()
	{
	// Choose a new rate
	u1 = rng->Uniform(FILE_AND_LINE);
	tmp_rate = u1*L;

	// Find which "rates" are on either side of the new rate. The quotes are there because the "rate" on the
	// left might be the lower bound (0.0) and the "rate" to the right might be the upper bound (L).
	tmp_rate_iter = std::lower_bound(model->gamma_rates_unnorm.begin(), model->gamma_rates_unnorm.end(), tmp_rate);
	if (tmp_rate_iter == model->gamma_rates_unnorm.begin())
		{
		rate_left = 0.0;
		rate_right = *tmp_rate_iter;
		}
	else if (tmp_rate_iter == model->gamma_rates_unnorm.end())
		{
		rate_left = *(tmp_rate_iter - 1);
		rate_right = L;
		}
	else
		{
		rate_left = *(tmp_rate_iter - 1);
		rate_right = *tmp_rate_iter;
		}

	// Locate probability corresponding to tmp_rate_iter
	unsigned offset = (unsigned)(tmp_rate_iter - model->gamma_rates_unnorm.begin());
	tmp_prob_iter = model->gamma_rate_probs.begin() + offset;

	// Choose a new category probability from the category probability prior distribution

	PHYCAS_ASSERT(cat_prob_prior);
	tmp_prob = cat_prob_prior->Sample();

	// Insert the new rate
	tmp_rate_iter = model->gamma_rates_unnorm.insert(tmp_rate_iter, tmp_rate);

	// Insert the new probability
	tmp_prob_iter = model->gamma_rate_probs.insert(tmp_prob_iter, tmp_prob);

	// Fix num_gamma_rates to reflect changes
	model->num_gamma_rates = (unsigned)model->gamma_rates_unnorm.size();

	PHYCAS_ASSERT(model->num_gamma_rates == model->gamma_rate_probs.size());
	PHYCAS_ASSERT(model->num_gamma_rates == ncat_after);

#if defined(SHOW_DEBUGGING_OUTPUT)
	std::cerr << "\n*********** In proposeAddCatMove() ***********" << std::endl;
	std::cerr << "  phi        = " << phi << std::endl;
	std::cerr << "  u1         = " << u1 << std::endl;
	std::cerr << "  tmp_rate   = " << tmp_rate << std::endl;
	std::cerr << "  tmp_prob   = " << tmp_prob << std::endl;
	std::cerr << "  rate_left  = " << rate_left << std::endl;
	std::cerr << "  rate_right = " << rate_right << std::endl;
	std::cerr << "  ncat_before = " << ncat_before << std::endl;
	std::cerr << "  ncat_after  = " << ncat_after << std::endl;
	std::cerr << "  model->num_gamma_rates         = " << model->num_gamma_rates << std::endl;
	std::cerr << "  model->gamma_rate_probs.size() = " << model->gamma_rate_probs.size() << std::endl;
	std::cerr << "  Unnormalized rates and probs after addcat:" << std::endl;
	for (unsigned i = 0; i < ncat_after; ++i)
		{
		std::cerr << str(boost::format("  %12.5f %12.5f") % model->gamma_rates_unnorm[i] % model->gamma_rate_probs[i])<< std::endl;
		}
	std::cerr << "~~~~~~~~~~~~~~~~~~" << std::endl;
#endif
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
void NCatMove::proposeDelCatMove()
	{
	// Select rate/prob to delete
	u1 = rng->Uniform(FILE_AND_LINE);
	unsigned which = (unsigned)(u1*ncat_before);

	// Save rate value, then remove from vector
	tmp_rate_iter = model->gamma_rates_unnorm.begin() + which;
	tmp_rate = *tmp_rate_iter;
	tmp_rate_iter = model->gamma_rates_unnorm.erase(tmp_rate_iter);
	if (tmp_rate_iter == model->gamma_rates_unnorm.begin())
		{
		rate_left = 0.0;
		rate_right = *tmp_rate_iter;
		}
	else if (tmp_rate_iter == model->gamma_rates_unnorm.end())
		{
		rate_left = *(tmp_rate_iter - 1);
		rate_right = L;
		}
	else
		{
		rate_left = *(tmp_rate_iter - 1);
		rate_right = *tmp_rate_iter;
		}

	// Save prob value, then remove from vector
	tmp_prob_iter = model->gamma_rate_probs.begin() + which;
	tmp_prob = *tmp_prob_iter;
	tmp_prob_iter = model->gamma_rate_probs.erase(tmp_prob_iter);

	// Fix num_gamma_rates to reflect changes
	model->num_gamma_rates = (unsigned)model->gamma_rates_unnorm.size();

	PHYCAS_ASSERT(model->num_gamma_rates == model->gamma_rate_probs.size());
	PHYCAS_ASSERT(model->num_gamma_rates == ncat_after);

#if defined(SHOW_DEBUGGING_OUTPUT)
	std::cerr << "\n*********** In proposeDelCatMove() ***********" << std::endl;
	std::cerr << "  phi        = " << phi << std::endl;
	std::cerr << "  u1         = " << u1 << std::endl;
	std::cerr << "  which      = " << which << std::endl;
	std::cerr << "  tmp_rate   = " << tmp_rate << std::endl;
	std::cerr << "  tmp_prob   = " << tmp_prob << std::endl;
	std::cerr << "  rate_left  = " << rate_left << std::endl;
	std::cerr << "  rate_right = " << rate_right << std::endl;
	std::cerr << "  Unnormalized rates and probs after delcat:" << std::endl;
	for (unsigned i = 0; i < ncat_after; ++i)
		{
		std::cerr << str(boost::format("  %12.5f %12.5f") % model->gamma_rates_unnorm[i] % model->gamma_rate_probs[i])<< std::endl;
		}
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Causes relative rates and probabilities to be returned to their state just 
|	prior to proposing the move.
*/
void NCatMove::revert()
	{
	MCMCUpdater::revert();
	if (addcat_move_proposed)
		{
		// Delete rate just inserted
		model->gamma_rates_unnorm.erase(tmp_rate_iter);

		// Delete prob just inserted
		model->gamma_rate_probs.erase(tmp_prob_iter);
		}
	else
		{
		// Reinsert rate at its original location
		model->gamma_rates_unnorm.insert(tmp_rate_iter, tmp_rate);

		// Reinsert prob at its original location
		model->gamma_rate_probs.insert(tmp_prob_iter, tmp_prob);
		}

	// Fix num_gamma_rates to reflect changes
	model->num_gamma_rates = (unsigned)model->gamma_rates_unnorm.size();

	PHYCAS_ASSERT(model->num_gamma_rates == model->gamma_rate_probs.size());
	PHYCAS_ASSERT(model->num_gamma_rates == ncat_before);

	// Renormalize the rates and probs so that we are ready for next likelihood calculation
	likelihood->recalcRelativeRates();

	reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
bool NCatMove::update()
	{
	bool accepted = false;
	total_updates++;

	if (is_fixed)
		return false;

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
	
#if DISABLED_UNTIL_WORKING_PRIOR_ACCOMMODATED
	double prev_ln_like = p->getLastLnLike();
	proposeNewState();

	likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
	double curr_ln_like		= likelihood->calcLnL(tree);
    double ln_like_ratio	= curr_ln_like - prev_ln_like;
    double ln_prior_ratio	= 0.0;

	if (addcat_move_proposed)
		{
		double ln_ncat_prior_ratio = std::log(lambda) - std::log((double)ncat_before);

		double logL = std::log(L);
		double float_s = (double)s;
		double s_plus_one = (double)s + 1.0;
		double nc = (double)ncat_before;
		double s_plus_one_times_nc_plus_two = s_plus_one*(nc + 2.0);
		double s_plus_one_times_nc_plus_one = s_plus_one*(nc + 1.0);

		double ln_rates_prior_ratio = 0.0;
		ln_rates_prior_ratio -= s_plus_one*logL;
		ln_rates_prior_ratio -= cdf.LnGamma(s_plus_one);
		ln_rates_prior_ratio += cdf.LnGamma(s_plus_one_times_nc_plus_two);
		ln_rates_prior_ratio -= cdf.LnGamma(s_plus_one_times_nc_plus_one);
		ln_rates_prior_ratio += float_s*std::log(tmp_rate - rate_left);
		ln_rates_prior_ratio += float_s*std::log(rate_right - tmp_rate);
		ln_rates_prior_ratio -= float_s*std::log(rate_right - rate_left);

		// Assuming new probs are drawn from their prior. In this case, the prior ratio is zero
		// because the density of the new prob in the numerator of the prob prior ratio is 
		// cancelled out by the density of the new prob that appears in the denominator of the
		// Jacobian
		double ln_probs_prior_ratio = 0.0;

		ln_prior_ratio = ln_ncat_prior_ratio + ln_rates_prior_ratio + ln_probs_prior_ratio;
		}
	else
		{
		double ln_ncat_prior_ratio = std::log((double)ncat_after) - std::log(lambda);

		double logL = std::log(L);
		double float_s = (double)s;
		double s_plus_one = (double)s + 1.0;
		double nc = (double)ncat_before;
		double s_plus_one_times_nc = s_plus_one*nc;
		double s_plus_one_times_nc_plus_one = s_plus_one*(nc + 1);

		double ln_rates_prior_ratio = 0.0;
		ln_rates_prior_ratio += s_plus_one*logL;
		ln_rates_prior_ratio += cdf.LnGamma(s_plus_one);
		ln_rates_prior_ratio += cdf.LnGamma(s_plus_one_times_nc);
		ln_rates_prior_ratio -= cdf.LnGamma(s_plus_one_times_nc_plus_one);
		ln_rates_prior_ratio += float_s*std::log(rate_right - rate_left);
		ln_rates_prior_ratio -= float_s*std::log(tmp_rate - rate_left);
		ln_rates_prior_ratio -= float_s*std::log(rate_right - tmp_rate);

		// Assuming new probs are drawn from their prior. In this case, the prior ratio is zero
		// because the density of the new prob in the denominator of the prob prior ratio is 
		// cancelled out by the density of the new prob that appears in the numerator of the
		// Jacobian
		double ln_probs_prior_ratio = 0.0;

		ln_prior_ratio = ln_ncat_prior_ratio + ln_rates_prior_ratio + ln_probs_prior_ratio;
		}

    ln_like_ratio *= heating_power;
    if (is_standard_heating)
        {
        ln_prior_ratio *= heating_power;
        }

	double ln_accept_ratio = ln_like_ratio + ln_prior_ratio + ln_hastings + ln_jacobian;
	double u = rng->Uniform(FILE_AND_LINE);

#if defined(SHOW_DEBUGGING_OUTPUT)
	std::cerr << "\n*********** In NCatMove::update() ***********" << std::endl;
	std::cerr << "  total_updates   = " << total_updates << std::endl;
	std::cerr << "  ln_like_ratio   = " << ln_like_ratio << std::endl;
	std::cerr << "  ln_prior_ratio  = " << ln_prior_ratio << std::endl;
	std::cerr << "  ln_hastings     = " << ln_hastings << std::endl;
	std::cerr << "  ln_jacobian     = " << ln_jacobian << std::endl;
	std::cerr << "  ln_accept_ratio = " << ln_accept_ratio << std::endl;
	std::cerr << "  u               = " << u << std::endl;
#endif

    accepted = (ln_accept_ratio >= 0.0 || std::log(u) <= ln_accept_ratio);

    if (save_debug_info)
        {
        debug_info = str(boost::format("NCatMove: rate %f, prob %f, %s, %s") % tmp_rate % tmp_prob % (addcat_move_proposed ? "addcat" : "delcat") % (accepted ? "accepted" : "rejected"));
        }

	if (accepted)
		{
#if defined(SHOW_DEBUGGING_OUTPUT)
		std::cerr << "  " << (addcat_move_proposed ? "ADDCAT" : "DELCAT") << " move accepted" << std::endl;
#endif

		p->setLastLnPrior(0.0);			//@POL do we really need to save the prior?
		p->setLastLnLike(curr_ln_like);
		accept();
		}
	else
		{
#if defined(SHOW_DEBUGGING_OUTPUT)
		std::cerr << "  " << (addcat_move_proposed ? "ADDCAT" : "DELCAT") << " move rejected" << std::endl;
#endif

		curr_ln_like	= p->getLastLnLike();
		curr_ln_prior	= p->getLastLnPrior();
		revert();
		}

#if defined(SHOW_DEBUGGING_OUTPUT)
	std::cerr << "  Unnormalized rates and probs just before leaving NCatMove::update():" << std::endl;
	unsigned curr_nrates = model->gamma_rates_unnorm.size();
	for (unsigned i = 0; i < curr_nrates; ++i)
		{
		std::cerr << str(boost::format("  %12.5f %12.5f") % model->gamma_rates_unnorm[i] % model->gamma_rate_probs[i])<< std::endl;
		}

	if (curr_nrates > ncat_distr.size())
		ncat_distr.resize(curr_nrates, 0);
	ncat_distr[curr_nrates - 1]++;
	std::cerr << "  Current ncat distribution:" << std::endl;
	unsigned ncat_distr_size = ncat_distr.size();
	unsigned total = std::accumulate(ncat_distr.begin(), ncat_distr.end(), 0, std::plus<unsigned>());
	for (unsigned j = 0; j < ncat_distr_size; ++j)
		{
		double proportion = (double)ncat_distr[j]/(double)total;
		std::cerr << str(boost::format("  %12d %12.5f") % (j+1) % proportion) << std::endl;
		}

	std::cerr << "\n" << std::endl;
#endif	//SHOW_DEBUGGING_OUTPUT
#endif	//DISABLED_UNTIL_WORKING_PRIOR_ACCOMMODATED
	return accepted;
	}
	
#endif	//FLEXCAT_MODEL

