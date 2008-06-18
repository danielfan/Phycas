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
#if POLPY_NEWWAY
#include "phycas/src/basic_tree_node.hpp"
#endif
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/unimap_nni_move.hpp"
#include "phycas/src/basic_tree.hpp"

#include "boost/format.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
bool UnimapNNIMove::update()
	{
	// The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
	// A proposed UnimapNNIMove cannot be accepted without changing edge lengths, so it is best to bail out now.
	if (is_fixed)
		return false;

    tree->renumberInternalNodes(tree->GetNTips()); //@POL this should be somewhere else

    ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
	double prev_ln_like = p->getLastLnLike();
	TreeNode * prev_likelihood_root = likelihood->getLikelihoodRoot();

	proposeNewState();

    //@POL this brute-force recalculation should be avoided for the sake of speed once the move is working
    likelihood->recalcSMatrix(tree);

    curr_ln_like = likelihood->calcLnL(tree);

	double prev_ln_prior = 0.0;
	if (star_tree_proposal)
		{
		prev_ln_prior		= p->calcExternalEdgeLenPriorUnnorm(orig_edge_len);

        double curr_edgelen = orig_node->GetEdgeLen();
		curr_ln_prior		= p->calcExternalEdgeLenPriorUnnorm(curr_edgelen);
		}
	else
		{
        PHYCAS_ASSERT(ndY->IsInternal());
        prev_ln_prior  = (ndX->IsInternal() ? p->calcInternalEdgeLenPriorUnnorm(x) : p->calcExternalEdgeLenPriorUnnorm(x));
        prev_ln_prior += p->calcInternalEdgeLenPriorUnnorm(y);
        prev_ln_prior += (ndZ->IsInternal() ? p->calcInternalEdgeLenPriorUnnorm(z) : p->calcExternalEdgeLenPriorUnnorm(z));

        double xnew = ndX->GetEdgeLen();
		curr_ln_prior  = (ndX->IsInternal() ? p->calcInternalEdgeLenPriorUnnorm(xnew) : p->calcExternalEdgeLenPriorUnnorm(xnew));

        double ynew = ndY->GetEdgeLen();
        curr_ln_prior += p->calcInternalEdgeLenPriorUnnorm(ynew);

        double znew = ndZ->GetEdgeLen();
        curr_ln_prior += (ndZ->IsInternal() ? p->calcInternalEdgeLenPriorUnnorm(znew) : p->calcExternalEdgeLenPriorUnnorm(znew));
		}

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

	double ln_accept_ratio = curr_posterior - prev_posterior + getLnHastingsRatio() + getLnJacobian();
    //double lnu = std::log(rng->Uniform(FILE_AND_LINE));
    //bool accepted = (ln_accept_ratio >= 0.0 || lnu <= ln_accept_ratio);
    //bool accepted = (ln_accept_ratio >= 0.0 || std::log(rng->Uniform(FILE_AND_LINE)) <= ln_accept_ratio);
    double lnu = DBL_MAX;
    bool accepted = (ln_accept_ratio >= 0.0);
    if (!accepted)
        {
        double u = rng->Uniform(FILE_AND_LINE);
        lnu = std::log(u);
        accepted = (lnu <= ln_accept_ratio);
        }

    if (save_debug_info)
        {
    	if (star_tree_proposal)
            {
            debug_info = str(boost::format("LS: %.5f -> %.5f (%s)") % orig_edge_len % orig_node->GetEdgeLen() % (accepted ? "accepted" : "rejected"));
            }
        else
            {
            debug_info = str(boost::format("UnimapNNIMove: topology %s, case = %d, x=%f, y=%f, z=%f, newX=%f, newY=%f, newZ=%f, %s, lnu = %.5f, lnr = %.5f, curr = %.5f, prev = %.5f") 
                % (topol_changed ? "changed" : "unchanged") 
                % which_case
                % x 
                % y 
                % z 
                % (ndX->GetEdgeLen()) 
                % (ndY->GetEdgeLen()) 
                % (ndZ->GetEdgeLen()) 
                % (accepted ? "accepted" : "rejected")
                % lnu
                % ln_accept_ratio
                % curr_posterior
                % prev_posterior);
            }
        }
    
    if (accepted)
		{
		p->setLastLnPrior(curr_ln_prior);
		p->setLastLnLike(curr_ln_like);
		accept();
		return true;
		}
	else
		{
		curr_ln_like	= p->getLastLnLike();
		curr_ln_prior	= p->getLastLnPrior();
		revert();

        //@POL 14-Mar-2008 First part of assert below added because prev_likelihood_root can legitimately be NULL
        // but it is troublesome that the original version of the assert was not tripped more often!
        // Why was it tripped only when running with no data? More thought should be given to this
        // situation.
		PHYCAS_ASSERT(!prev_likelihood_root || prev_likelihood_root->IsInternal());
		likelihood->useAsLikelihoodRoot(prev_likelihood_root);
		return false;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapNNIMove::proposeNewState()
    {
    }

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapNNIMove::revert()
	{
    }

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void UnimapNNIMove::accept()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The default constructor sets `lambda' to the default value (0.2), sets `topol_changed' to false, and `m' and `mstar'
|	to 0.0. All other data members are automatically initialized (shared pointers) or are initialized via a call to 
|	reset().
*/
UnimapNNIMove::UnimapNNIMove() : MCMCUpdater()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Forgets information saved to enable reverting a proposed move.
*/
void UnimapNNIMove::reset()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move.
*/
double UnimapNNIMove::getLnHastingsRatio() const
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Jacobian for this move.
*/
double UnimapNNIMove::getLnJacobian() const
	{
	return 0.0;
	}

}	// namespace phycas

#endif
