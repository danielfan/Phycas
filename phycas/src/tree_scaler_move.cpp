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
#if POLPY_NEWWAY
	is_move = true;
    n = 0;
    m = 0.0;
    mstar = 0.0;
    forward_scaler = 1.0;
    reverse_scaler = 1.0;
    lambda = 0.5;
#else
	curr_value = 1.0;
	has_slice_sampler = true;
	is_move = true;
	is_master_param = false;
	is_hyper_param = false;
#endif
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'lambda', which is the tuning parameter for this move.
*/
void TreeScalerMove::setLambda(
  double x) /* is the new value for `lambda' */
	{
	lambda = x;
	}
#endif

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the data member 'lambda', which is the tuning parameter for this move.
*/
double TreeScalerMove::getLambda() const
	{
	return lambda;
	}
#endif

#if POLPY_NEWWAY
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
#endif

#if POLPY_NEWWAY
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
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool TreeScalerMove::update()
	{
    if (is_fixed)
		return false;

#if POLPY_NEWWAY
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

#else

    //std::cerr << "~~~~~~ TreeScalerMove::update()" << std::endl;

    slice_sampler->Sample();

    //double v = slice_sampler->GetLastSampledXValue();
    //std::ofstream outf("treescaler.txt", std::ios::out | std::ios::app);
    //outf.setf(std::ios::floatfield, std::ios::fixed);
    //outf.setf(std::ios::showpoint);
    //outf << v << std::endl;
    //outf.close();

    if (save_debug_info)
        {
        debug_info = str(boost::format("TreeScalerMove %f") % (slice_sampler->GetLastSampledXValue()));
        }

    rescaleAllEdgeLengths();
#endif

    return true;
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Reverses move made in proposeNewState.
*/
void TreeScalerMove::revert()
	{
    tree->ScaleAllEdgeLens(reverse_scaler);
	}
#endif

#if POLPY_NEWWAY
/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void TreeScalerMove::accept()
	{
	}
#endif

#if POLPY_NEWWAY
#else
/*----------------------------------------------------------------------------------------------------------------------
|	TreeScalerMove is a functor whose operator() returns a value proportional to the full-conditional posterior  
|	probability density for a particular scaling of the entire tree. If the supplied scaling factor (f) is out of bounds
|   (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative infinity).
*/
double TreeScalerMove::operator()(
  double f)	/**< is a new value for the scaling factor */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	if (f > 0.0)
		{
		curr_value = f;
        tree->SetTreeScale(f);
		recalcPrior();

        likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
		curr_ln_like = likelihood->calcLnL(tree);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

        if (is_standard_heating)
            return heating_power*(curr_ln_like + curr_ln_prior);
        else
            return heating_power*curr_ln_like + curr_ln_prior;
		}
    else
        return ln_zero;
	}
#endif

#if POLPY_NEWWAY
#else
/*----------------------------------------------------------------------------------------------------------------------
|	Scales all edge lengths in the tree so that the tree is now different from its original length by a factor of 
|   `curr_value'.
*/
void TreeScalerMove::rescaleAllEdgeLengths()
	{
    //@POL this function is replicated in TreeManip - should just use that version
    double scaling_factor = curr_value;

    curr_value = 1.0;
    tree->SetTreeScale(1.0);
    slice_sampler->SetXValue(1.0);

    // Change the edge lengths
	for (TreeNode * nd = tree->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		// All nodes have an edge associated with them except for the root
		//
		if (!nd->IsTipRoot())
			{
            double old_edgelen = nd->GetEdgeLen();
            double new_edgelen = old_edgelen*scaling_factor;
            nd->SetEdgeLen(new_edgelen);
			}
		}
	}
#endif

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
