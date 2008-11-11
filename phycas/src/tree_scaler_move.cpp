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
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'lambda', which is the tuning parameter for this move.
*/
void TreeScalerMove::setLambda(
  double x) /* is the new value for `lambda' */
	{
	lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the data member 'lambda', which is the tuning parameter for this move.
*/
double TreeScalerMove::getLambda() const
	{
	return lambda;
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

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);

    double prev_ln_prior		= recalcPrior();
	double prev_ln_like			= p->getLastLnLike();

    proposeNewState();

    likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
	double curr_ln_prior		= recalcPrior();
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

    return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reverses move made in proposeNewState.
*/
void TreeScalerMove::revert()
	{
    tree->ScaleAllEdgeLens(reverse_scaler);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void TreeScalerMove::accept()
	{
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

}	// namespace phycas
