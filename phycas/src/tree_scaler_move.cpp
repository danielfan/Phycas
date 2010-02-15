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
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/tree_scaler_move.hpp"
#include "phycas/src/basic_tree.hpp"

#include "boost/format.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	The default constructor.
*/
TreeScalerMove::TreeScalerMove() : MCMCUpdater()
	{
	is_move = true;
    n = 0;
    m = 0.0;
    mstar = 0.0;
    forward_scaler = 1.0;
    reverse_scaler = 1.0;
    lambda = 0.5;
    min_lambda = 0.5;
    max_lambda = 0.5;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'min_lambda', which is the tuning parameter used for exploring the posterior 
|   distribution in this move.
*/
void TreeScalerMove::setPosteriorTuningParam(
  double x) /* is the new value for `min_lambda' */
	{
	min_lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member `max_lambda', which is the tuning parameter used for exploring the prior 
|   distribution in this move.
*/
void TreeScalerMove::setPriorTuningParam(
  double x) /* is the new value for `max_lambda' */
	{
	max_lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'lambda', which is the tuning parameter for this move, based on a boldness value
|	that ranges from 0 (least bold) to 100 (most bold). Simple linear interpolation is used (i.e. a boldness of 50
|   results in `lambda' halfway between `min_lambda' and `max_lambda').
*/
void TreeScalerMove::setBoldness(
  double x) /* is the new boldness value */
	{
	boldness = x;
	if (boldness < 0.0)
		boldness = 0.0;
	else if (boldness > 100.0)
		boldness = 100.0;

    // compute lambda from boldness value
	lambda = min_lambda + (max_lambda - min_lambda)*boldness/100.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move. The Hastings ratio is (`mstar'/`m')^n, where `mstar' is
|	the new tree length and `m' is the tree length before the move is proposed. The value n is the number of edges in
|   the tree.
*/
double TreeScalerMove::getLnHastingsRatio() const
	{
    double ln_hastings = (double)n*(std::log(mstar) - std::log(m));
    return ln_hastings;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Proposes a new value by which to scale the tree length.
*/
void TreeScalerMove::proposeNewState()
	{
    n = tree->GetNNodes() - 1;  // root node's edge does not count
    m = tree->EdgeLenSum();
    mstar = m*exp(lambda*(rng->Uniform(FILE_AND_LINE) - 0.5));
    forward_scaler = mstar/m;
    reverse_scaler = m/mstar;
    tree->ScaleAllEdgeLens(forward_scaler);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool TreeScalerMove::update()
	{
    if (is_fixed)
		return false;
		
	// std::cerr << "****** TreeScalerMove::update" << std::endl;

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);

    double prev_ln_prior			= recalcPrior();
	double prev_ln_like				= p->getLastLnLike();
#if POLPY_NEWWAY
	double prev_ln_working_prior	= (use_working_prior ? recalcWorkingPriorForMove() : 0.0);
#endif

    proposeNewState();

    likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
	double curr_ln_prior			= recalcPrior();
	double curr_ln_like				= (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
#if POLPY_NEWWAY
	double curr_ln_working_prior	= (use_working_prior ? recalcWorkingPriorForMove() : 0.0);
#endif

    double prev_posterior = 0.0;
	double curr_posterior = 0.0;
    if (is_standard_heating)
        {
        prev_posterior = heating_power*(prev_ln_like + prev_ln_prior);
	    curr_posterior = heating_power*(curr_ln_like + curr_ln_prior);
#if POLPY_NEWWAY
		if (use_working_prior)
			{
			prev_posterior += (1.0 - heating_power)*prev_ln_working_prior;
			curr_posterior += (1.0 - heating_power)*curr_ln_working_prior;
			}
#endif
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
	    if (save_debug_info)
    	    {
			debug_info = boost::str(boost::format("ACCEPT, forward_scaler = %.5f, prev_ln_prior = %.5f, curr_ln_prior = %.5f, prev_ln_like = %.5f, curr_ln_like = %.5f, lnu = %.5f, ln_accept_ratio = %.5f") % forward_scaler % prev_ln_prior % curr_ln_prior % prev_ln_like % curr_ln_like % lnu % ln_accept_ratio);
			}
		p->setLastLnPrior(curr_ln_prior);
		p->setLastLnLike(curr_ln_like);
		accept();
		return true;
		}
	else
		{
	    if (save_debug_info)
    	    {
			debug_info = boost::str(boost::format("REJECT, forward_scaler = %.5f, prev_ln_prior = %.5f, curr_ln_prior = %.5f, prev_ln_like = %.5f, curr_ln_like = %.5f, lnu = %.5f, ln_accept_ratio = %.5f") % forward_scaler % prev_ln_prior % curr_ln_prior % prev_ln_like % curr_ln_like % lnu % ln_accept_ratio);
			}
		curr_ln_like	= p->getLastLnLike();
		curr_ln_prior	= p->getLastLnPrior();
		revert();
		return false;
		}

    return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reverses move made in proposeNewState.
*/
void TreeScalerMove::revert()
	{
	MCMCUpdater::revert();
    tree->ScaleAllEdgeLens(reverse_scaler);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void TreeScalerMove::accept()
	{
	MCMCUpdater::accept();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the joint log prior over all edges in the associated tree and sets `curr_ln_prior'.
*/
double TreeScalerMove::recalcPrior()
	{	
    // Loop through all EdgeLenMasterParam objects and call the recalcPrior function of each.
    // Each EdgeLenMasterParam object knows how to compute the prior for the edge lengths it controls.
    curr_ln_prior = 0.0;
	ChainManagerShPtr p = chain_mgr.lock();
    const MCMCUpdaterVect & edge_length_params = p->getEdgeLenParams();
    for (MCMCUpdaterVect::const_iterator it = edge_length_params.begin(); it != edge_length_params.end(); ++it)
        {
        curr_ln_prior += (*it)->recalcPrior();
        }

	return curr_ln_prior;
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Computes the joint log working prior over all edges in the associated tree.
*/
double TreeScalerMove::recalcWorkingPriorForMove() const
	{	
    // Loop through all EdgeLenMasterParam objects and call the recalcWorkingPrior function of each.
    // Each EdgeLenMasterParam object knows how to compute the working prior for the edge lengths it controls.
    double ln_working_prior = 0.0;
	ChainManagerShPtr p = chain_mgr.lock();
    const MCMCUpdaterVect & edge_length_params = p->getEdgeLenParams();
    for (MCMCUpdaterVect::const_iterator it = edge_length_params.begin(); it != edge_length_params.end(); ++it)
        {
        ln_working_prior += (*it)->recalcWorkingPrior();
        }

	return ln_working_prior;
	}
#endif

}	// namespace phycas
