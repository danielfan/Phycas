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
#include "phycas/src/edge_move.hpp"
#include "phycas/src/topo_prior_calculator.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/tree_manip.hpp"

//#define DEBUG_LOG
#if defined(DEBUG_LOG)
#   include <fstream>
#endif

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor sets `origEdgelen' to 0.0, `origNode' to NULL, `likeRoot' to NULL, and `lambda' to 1.0.
*/
EdgeMove::EdgeMove() : MCMCUpdater()
	{
	is_move = true;
	lambda	= 1.0;
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Sets `origNode' to NULL and `origEdgelen' to 0.0. These variables are used to save quantities required for reverting
|	a proposed move that was not accepted, so this function is called at the end of both accept() and revert() to reset the
|	object for the next proposal. It is also called by the constructor to initialize these variables.
*/
void EdgeMove::reset()
	{
	origEdgelen	= 0.0;
	origNode	= NULL;
	likeRoot	= NULL;

	// one_edgelen should have one element, used for computing the edge length prior in update
	if (one_edgelen.empty())
		one_edgelen.push_back(0.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'min_lambda', which is the tuning parameter used for exploring the posterior 
|   distribution in this move.
*/
void EdgeMove::setPosteriorTuningParam(
  double x) /* is the new value for `min_lambda' */
	{
	min_lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member `max_lambda', which is the tuning parameter used for exploring the prior 
|   distribution in this move.
*/
void EdgeMove::setPriorTuningParam(
  double x) /* is the new value for `max_lambda' */
	{
	max_lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'lambda', which is the tuning parameter for this move, based on a boldness value
|	that ranges from 0 (least bold) to 100 (most bold). Simple linear interpolation is used (i.e. a boldness of 50
|   results in `lambda' halfway between `min_lambda' and `max_lambda').
*/
void EdgeMove::setBoldness(
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
|	Sets the value for the data member lambda, which is the tuning parameter for this move.
*/
void EdgeMove::setLambda(double x)
	{
	lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the data member lambda, which is the tuning parameter for this move.
*/
double EdgeMove::getLambda() const
	{
	return lambda;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move. The Hastings ratio is:
|>
|	Pr(new state -> old state)   [1/(lambda m)]
|	-------------------------- = --------------- = m* / m
|	Pr(old state -> new state)   [1/(lambda m*)]
|>
|	Here, m* is the new length of the modified edge and m is the length before the move is proposed.
*/
double EdgeMove::getLnHastingsRatio() const
	{
	return std::log(origNode->GetEdgeLen()) - std::log(origEdgelen);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	This move does not change the model dimension, so the Jacobian is irrelevant.
*/
double EdgeMove::getLnJacobian() const
	{
	return 0.0;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is accepted.
*/
void EdgeMove::accept()
	{
	MCMCUpdater::accept();
    PHYCAS_ASSERT(likeRoot);
    likelihood->useAsLikelihoodRoot(likeRoot);
    likelihood->invalidateAwayFromNode(*likeRoot);
    likelihood->invalidateBothEnds(origNode);	//@POL really just need invalidateParentalOnly function

	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Causes tree to be returned to its state just prior to proposing the move.
*/
void EdgeMove::revert()
	{
	MCMCUpdater::revert();
	origNode->SetEdgeLen(origEdgelen);

    PHYCAS_ASSERT(likeRoot);
    likelihood->useAsLikelihoodRoot(likeRoot);
    likelihood->invalidateAwayFromNode(*likeRoot);
    likelihood->invalidateBothEnds(origNode);	//@POL really just need invalidateParentalOnly function

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
	for (origNode = tree->GetFirstPreorder(); origNode != NULL; origNode = origNode->GetNextPreorder())
		{
		// All nodes have an edge associated with them except for the root
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

    //origNode->SelectNode();
    //likelihood->startTreeViewer(tree, "EdgeMove selected edge");
    //origNode->UnselectNode();

	// Modify the edge
	double m		= origNode->GetEdgeLen();
	double mstar	= m*std::exp(lambda*(rng->Uniform(FILE_AND_LINE) - 0.5));
	origNode->SetEdgeLen(mstar);

    if (origNode->IsInternal())
        likeRoot = origNode;
    else
        likeRoot = origNode->GetParent();
    PHYCAS_ASSERT(likeRoot);
    likelihood->useAsLikelihoodRoot(likeRoot);
    likelihood->invalidateAwayFromNode(*likeRoot);
    likelihood->invalidateBothEnds(origNode);	//@POL really just need invalidateParentalOnly function
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
bool EdgeMove::update()
	{
	//need to SAMCize
	
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
	double prev_ln_working_prior = 0.0;
#if USING_EDGE_SPECIFIC_WORKING_PRIORS
	if (use_working_prior)
		prev_ln_working_prior = (is_internal_edge ? p->calcInternalEdgeLenWorkingPrior(*origNode, origEdgelen) : p->calcExternalEdgeLenWorkingPrior(*origNode, origEdgelen));
#else
	if (use_working_prior)
		prev_ln_working_prior = (is_internal_edge ? p->calcInternalEdgeLenWorkingPrior(origEdgelen) : p->calcExternalEdgeLenWorkingPrior(origEdgelen));
#endif

	double curr_ln_like			= (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);

	//one_edgelen[0]				= origNode->GetEdgeLen();
	//double curr_ln_prior		= p->partialEdgeLenPrior(one_edgelen);
    double curr_edgelen         = origNode->GetEdgeLen();
	double curr_ln_prior		= (is_internal_edge ? p->calcInternalEdgeLenPriorUnnorm(curr_edgelen) : p->calcExternalEdgeLenPriorUnnorm(curr_edgelen));
	double curr_ln_working_prior = 0.0;
#if USING_EDGE_SPECIFIC_WORKING_PRIORS
	if (use_working_prior)
		curr_ln_working_prior = (is_internal_edge ? p->calcInternalEdgeLenWorkingPrior(*origNode, curr_edgelen) : p->calcExternalEdgeLenWorkingPrior(*origNode, curr_edgelen));
#else
	if (use_working_prior)
		curr_ln_working_prior = (is_internal_edge ? p->calcInternalEdgeLenWorkingPrior(curr_edgelen) : p->calcExternalEdgeLenWorkingPrior(curr_edgelen));
#endif

    double prev_posterior = 0.0;
	double curr_posterior = 0.0;
#if defined(DEBUG_LOG)
    std::ofstream doof("doof.txt", std::ios::app);
#endif
    if (is_standard_heating)
        {
        prev_posterior = heating_power*(prev_ln_like + prev_ln_prior);
	    curr_posterior = heating_power*(curr_ln_like + curr_ln_prior);
		if (use_working_prior)
			{
			prev_posterior += (1.0 - heating_power)*prev_ln_working_prior;
			curr_posterior += (1.0 - heating_power)*curr_ln_working_prior;
			}
#if defined(DEBUG_LOG)
        doof << boost::str(boost::format("std\t%.5f\t%.5f\t%.5f\t%.5f\t") % heating_power % curr_ln_like % curr_ln_prior % curr_posterior);
#endif
        }
    else
        {
        prev_posterior = heating_power*prev_ln_like + prev_ln_prior;
	    curr_posterior = heating_power*curr_ln_like + curr_ln_prior;
#if defined(DEBUG_LOG)
        doof << boost::str(boost::format("like\t%.5f\t%.5f\t%.5f\t%.5f\t") % heating_power % curr_ln_like % curr_ln_prior % curr_posterior);
#endif
        }

	double ln_hastings			= getLnHastingsRatio();

	double ln_accept_ratio		= curr_posterior - prev_posterior + ln_hastings;

    double lnu = std::log(rng->Uniform(FILE_AND_LINE));
	if (ln_accept_ratio >= 0.0 || lnu <= ln_accept_ratio)
		{
		p->setLastLnPrior(curr_ln_prior);
		p->setLastLnLike(curr_ln_like);
#if defined(DEBUG_LOG)
        doof << boost::str(boost::format("%.5f\t%.5f\t%.5f\t%.5f\t%s") % origEdgelen % curr_edgelen % lnu % ln_accept_ratio % "accept") << std::endl;
#endif
		accept();
		return true;
		}
	else
		{
		curr_ln_like	= p->getLastLnLike();
		curr_ln_prior	= p->getLastLnPrior();
#if defined(DEBUG_LOG)
        doof << boost::str(boost::format("%.5f\t%.5f\t%.5f\t%.5f\t%s") % origEdgelen % curr_edgelen % lnu % ln_accept_ratio % "reject") << std::endl;
#endif
		revert();
		return false;
		}
#if defined(DEBUG_LOG)
    doof.close();
#endif
	}

