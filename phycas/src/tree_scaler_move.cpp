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

#include "phycas/src/cipres/CipresDataMatrixHelper.h"
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
	curr_value = 1.0;
    prev_value = 1.0;
	has_slice_sampler = true;
	is_move = true;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
void TreeScalerMove::update()
	{
	if (is_fixed)
		return;

    // Each time the update function is called, we assume that the current scale is 1.0
    // The Sample() function will call operator() many times with new curr_value values
    // each of which will be relative to the starting scale.
    curr_value = 1.0;
    prev_value = 1.0;
    slice_sampler->Sample();

    if (save_debug_info)
        {
        debug_info = str(boost::format("TreeScalerMove %f") % (slice_sampler->GetLastSampledXValue()));
        }
	}

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
        prev_value = curr_value;
		curr_value = f;
        scaleAllEdgeLengths();
		recalcPrior();

		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
		curr_ln_like = likelihood->calcLnL(tree);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);
		}
	return curr_ln_like + curr_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Scales all edge lengths in the tree so that the tree is now different from its original length by a factor of 
|   `curr_value'.
*/
void TreeScalerMove::scaleAllEdgeLengths()
	{
    // Suppose we have a tree length (TL) equal to 10 at the start.
    // In the first call to operator(), the tree is scaled by factor 0.5, yielding TL = 5
    // In the next call to operator(), suppose the scaling factor (curr_value) becomes 0.2
    // We cannot simply multiply all edges in the tree by 0.2 because the tree has been scaled previously
    // and is now half as long as it was at the start. Need to scale by 0.2/0.5 (curr_value/prev_value),
    // which yields (5)(0.2)/0.5 = 2.
    double scaling_factor = curr_value/prev_value;

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
