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

#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/basic_tree_node.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/state_freq_move.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/tree_manip.hpp"

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor sets tuning parameter `psi' to 300.0 and calls reset().
*/
StateFreqMove::StateFreqMove() : MCMCUpdater()
	{
	psi = 300.0;
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|   Resets the 'dir_forward' and 'dir_reverse' shared pointers.
*/
void StateFreqMove::reset()
	{
    dir_forward.reset();
    dir_reverse.reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'psi', which is the tuning parameter for this move.
*/
void StateFreqMove::setPsi(
  double x) /* is the new value for psi */
	{
	psi = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the data member 'psi', which is the tuning parameter for this move.
*/
double StateFreqMove::getPsi() const
	{
	return psi;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|   Chooses a new vector of state frequencies using a sharp Dirichlet distribution centered at the original frequencies.
*/
void StateFreqMove::proposeNewState()
	{
    // copy the current frequencies to the data member orig_freqs
    const std::vector<double> & rfreqs = model->getStateFreqs();
    orig_freqs.resize(rfreqs.size());
    std::copy(rfreqs.begin(), rfreqs.end(), orig_freqs.begin());

    // create vector of Dirichlet parameters for selecting new frequencies
    c_forward.resize(rfreqs.size());
	std::transform(orig_freqs.begin(), orig_freqs.end(), c_forward.begin(), boost::lambda::_1*psi);

    // create Dirichlet distribution and sample from it to get proposed frequencies
    dir_forward = DirichletShPtr(new DirichletDistribution(c_forward));
    dir_forward->SetLot(getLot().get());
    new_freqs = dir_forward->Sample();

    // create vector of Dirichlet parameters for selecting old frequencies (needed for Hasting ratio calculation)
    c_reverse.resize(new_freqs.size());
	std::transform(new_freqs.begin(), new_freqs.end(), c_reverse.begin(), boost::lambda::_1*psi);
    dir_reverse = DirichletShPtr(new DirichletDistribution(c_reverse));

    //for (unsigned i = 0; i < 5; ++i)
    //    {
    //    std::vector<double> tmp_freqs = dir_forward->Sample();
    //    std::cerr << boost::str(boost::format("  %d %.9f %.9f %.9f %.9f\n") % tmp_freqs[0] % i % tmp_freqs[1] % tmp_freqs[2] % tmp_freqs[3]);
    //    }
    //std::cerr << boost::str(boost::format("%.9f %.9f %.9f %.9f ") % new_freqs[0] % new_freqs[1] % new_freqs[2] % new_freqs[3]);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move. The Hastings ratio is:
|>
|	Pr(new state -> old state)   Dir(pi'[i]*psi) density evaluated at pi
|	-------------------------- = ---------------------------------------
|	Pr(old state -> new state)   Dir(pi[i]*psi) density evaluated at pi'
|>
|   where pi[i] is the ith current frequency and pi'[i] is the ith proposed new frequency.
*/
double StateFreqMove::getLnHastingsRatio() const
	{
    return (dir_reverse->GetLnPDF(orig_freqs) - dir_forward->GetLnPDF(new_freqs));
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	This move does not change the model dimension, so the Jacobian is irrelevant.
*/
double StateFreqMove::getLnJacobian() const
	{
	return 0.0;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is accepted. Simply calls the reset() function.
*/
void StateFreqMove::accept()
	{
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Calls the reset() function after reinstating the original state frequencies
|   and ensuring that all conditional likelihood arrays will be recalculated when the likelihood is next calculated.
*/
void StateFreqMove::revert()
	{
    // replace current frequencies with original ones
    model->setStateFreqsUnnorm(orig_freqs);

    // invalidate all CLAs
    likelihood->useAsLikelihoodRoot(NULL);	

	reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
bool StateFreqMove::update()
	{
	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);

	proposeNewState();

    double prev_ln_prior		= mvprior->GetRelativeLnPDF(orig_freqs);
	double prev_ln_like			= p->getLastLnLike();

    // replace current frequencies with new ones
    model->setStateFreqsUnnorm(new_freqs);
    likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs

	double curr_ln_prior		= mvprior->GetRelativeLnPDF(new_freqs);
	double curr_ln_like			= likelihood->calcLnL(tree);

    double prev_posterior = 0.0;
	double curr_posterior = 0.0;
    if (is_standard_heating)
        {
        prev_posterior = heating_power*(prev_ln_like + prev_ln_prior);
	    curr_posterior = heating_power*(curr_ln_like + curr_ln_prior);
        }
    else
        {
        prev_posterior = heating_power*prev_ln_like + prev_ln_prior;
	    curr_posterior = heating_power*curr_ln_like + curr_ln_prior;
        }

	double ln_hastings			= getLnHastingsRatio();
	double ln_accept_ratio		= curr_posterior - prev_posterior + ln_hastings;

    double lnu = std::log(rng->Uniform(FILE_AND_LINE));
	if (ln_accept_ratio >= 0.0 || lnu <= ln_accept_ratio)
		{
        //std::cerr << "ACCEPTED\n" << std::endl;
		p->setLastLnPrior(curr_ln_prior);
		p->setLastLnLike(curr_ln_like);
		accept();
		return true;
		}
	else
		{
        //std::cerr << "rejected\n" << std::endl;
		curr_ln_like	= p->getLastLnLike();
		curr_ln_prior	= p->getLastLnPrior();
		revert();
		return false;
		}
	}
