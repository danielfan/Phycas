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
#include "phycas/src/cipres/CipresDataMatrixHelper.h"
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/edge_move.hpp"
#include "phycas/src/topo_prior_calculator.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/tree_manip.hpp"

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor sets `origEdgelen' to 0.0, `origNode' to NULL, and `lambda' to 1.0.
*/
EdgeMove::EdgeMove() : MCMCUpdater()
	{
	lambda			= 1.0;
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is accepted.
*/
void EdgeMove::accept()
	{
	likelihood->useAsLikelihoodRoot(origNode);
	likelihood->discardCacheAwayFromNode(*origNode);
	likelihood->discardCacheBothEnds(origNode);

	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Causes tree to be returned to its state just prior to proposing the move.
*/
void EdgeMove::revert()
	{
	origNode->SetEdgeLen(origEdgelen);

	likelihood->useAsLikelihoodRoot(origNode);
	likelihood->restoreFromCacheAwayFromNode(*origNode);
	likelihood->restoreFromCacheParentalOnly(origNode);

	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Chooses a random edge and changes its current length m to a new length m* using the following formula, where `lambda' is
|	a tuning parameter.
|>
|	m* = m*exp(lambda*(r.Uniform(FILE_AND_LINE) - 0.5))
|>
*/
void EdgeMove::proposeNewState()
	{
	// Choose edge randomly.
	//
	unsigned numEdges = tree->GetNNodes() - 1;
	unsigned k = rng->SampleUInt(numEdges);
	unsigned i = 0;
	//@POL this loop is crying out for the for_each algorithm
	for (origNode = tree->GetFirstPreorder(); origNode != NULL; origNode = origNode->GetNextPreorder())
		{
		// All nodes have an edge associated with them except for the root
		//
		if (!origNode->IsTipRoot())
			{
			if (i == k)
				{
				origEdgelen = origNode->GetEdgeLen();
				break;
				}
			++i;
			}
		}

	// Modify the edge
	//
	double m		= origNode->GetEdgeLen();
	double mstar	= m*std::exp(lambda*(rng->Uniform(FILE_AND_LINE) - 0.5));
	origNode->SetEdgeLen(mstar);

	likelihood->useAsLikelihoodRoot(origNode);
	likelihood->invalidateAwayFromNode(*origNode);
	likelihood->invalidateBothEnds(origNode);	//@POL really just need invalidateParentalOnly function
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
bool EdgeMove::update()
	{
	// The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
	// A proposed EdgeMove cannot be accepted without changing edge lengths, so it is best to just bail out now.
	if (is_fixed)
		return false;

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);

	double prev_ln_like			= p->getLastLnLike();

	proposeNewState();

	//one_edgelen[0]				= origEdgelen;
	//double prev_ln_prior		= p->partialEdgeLenPrior(one_edgelen);
    bool is_internal_edge       = origNode->IsInternal();
    double prev_ln_prior		= (is_internal_edge ? p->calcInternalEdgeLenPriorUnnorm(origEdgelen) : p->calcExternalEdgeLenPriorUnnorm(origEdgelen));

	double curr_ln_like			= likelihood->calcLnL(tree);

	//one_edgelen[0]				= origNode->GetEdgeLen();
	//double curr_ln_prior		= p->partialEdgeLenPrior(one_edgelen);
    double curr_edgelen         = origNode->GetEdgeLen();
	double curr_ln_prior		= (is_internal_edge ? p->calcInternalEdgeLenPriorUnnorm(curr_edgelen) : p->calcExternalEdgeLenPriorUnnorm(curr_edgelen));

#if POLPY_NEWWAY
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
#else
	double prev_posterior = prev_ln_like + prev_ln_prior;
	double curr_posterior = curr_ln_like + curr_ln_prior;
#endif

	double ln_hastings			= getLnHastingsRatio();

	double ln_accept_ratio		= curr_posterior - prev_posterior + ln_hastings;

	if (ln_accept_ratio >= 0.0 || std::log(rng->Uniform(FILE_AND_LINE)) <= ln_accept_ratio)
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
		return false;
		}
	}

