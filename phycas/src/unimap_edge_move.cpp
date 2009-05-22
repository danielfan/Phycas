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
#include "phycas/src/unimap_edge_move.hpp"
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
UnimapEdgeMove::UnimapEdgeMove() : MCMCUpdater()
	{
	lambda			= 1.0;
	mdot	= 0;
	r	= 1.0;
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Sets `origNode' to NULL and `origEdgelen' to 0.0. These variables are used to save quantities required for reverting
|	a proposed move that was not accepted, so this function is called at the end of both accept() and revert() to reset the
|	object for the next proposal. It is also called by the constructor to initialize these variables.
*/
void UnimapEdgeMove::reset()
	{
	origEdgelen	= 0.0;
	origNode	= NULL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'min_lambda', which is the tuning parameter used for exploring the posterior 
|   distribution in this move.
*/
void UnimapEdgeMove::setPosteriorTuningParam(
  double x) /* is the new value for `min_lambda' */
	{
	min_lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member `max_lambda', which is the tuning parameter used for exploring the prior 
|   distribution in this move.
*/
void UnimapEdgeMove::setPriorTuningParam(
  double x) /* is the new value for `max_lambda' */
	{
	max_lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'lambda', which is the tuning parameter for this move, based on a boldness value
|	that ranges from 0 (least bold) to 100 (most bold). Simple linear interpolation is used (i.e. a boldness of 50
|   results in `lambda' halfway between `min_lambda' and `max_lambda').
*/
void UnimapEdgeMove::setBoldness(
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
void UnimapEdgeMove::setLambda(double x)
	{
	lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the data member lambda, which is the tuning parameter for this move.
*/
double UnimapEdgeMove::getLambda() const
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
double UnimapEdgeMove::getLnHastingsRatio() const
	{
	return std::log(r);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	This move does not change the model dimension, so the Jacobian is irrelevant.
*/
double UnimapEdgeMove::getLnJacobian() const
	{
	return 0.0;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is accepted.
*/
void UnimapEdgeMove::accept()
	{
	origNode->SetEdgeLen(origEdgelen*r);
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Causes tree to be returned to its state just prior to proposing the move.
*/
void UnimapEdgeMove::revert()
	{
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Chooses a random edge and changes its current length m to a new length m* using the following formula, where `lambda' is
|	a tuning parameter.
|>
|	m* = m*exp(lambda*(r.Uniform(FILE_AND_LINE) - 0.5))
|>
*/
void UnimapEdgeMove::proposeNewState()
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

	// Modify the edge
	const Univents & u =  getUniventsConstRef(*origNode);
	if (!u.isValid())
		{
		likelihood->GetUniventProbMgrRef().recalcUMat();
		likelihood->remapUniventsForNode(tree, origNode);
		}
	mdot = u.getMDot();
	r = std::exp(lambda*(rng->Uniform(FILE_AND_LINE) - 0.5));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
bool UnimapEdgeMove::update()
	{
	// The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
	// A proposed UnimapEdgeMove cannot be accepted without changing edge lengths, so it is best to just bail out now.
	if (is_fixed)
		return false;

//	std::cerr << "****** UnimapEdgeMove::update" << std::endl;

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);

	proposeNewState();
	double nsites = (double)likelihood->getNPatterns();
	double uniformization_lambda = model->calcUniformizationLambda();
	
    bool is_internal_edge       = origNode->IsInternal();
    double prev_ln_prior		= (is_internal_edge ? p->calcInternalEdgeLenPriorUnnorm(origEdgelen) : p->calcExternalEdgeLenPriorUnnorm(origEdgelen));

    double curr_edgelen         = r*origEdgelen;
	double curr_ln_prior		= (is_internal_edge ? p->calcInternalEdgeLenPriorUnnorm(curr_edgelen) : p->calcExternalEdgeLenPriorUnnorm(curr_edgelen));

	double log_posterior_ratio = 0.0;
	const double log_prior_ratio = curr_ln_prior - prev_ln_prior;
    const double log_likelihood_ratio = (double)mdot*log(r) - nsites*uniformization_lambda*(curr_edgelen - origEdgelen);
    if (is_standard_heating)
        log_posterior_ratio = heating_power*(log_likelihood_ratio + log_prior_ratio);
    else
	    log_posterior_ratio = heating_power*log_likelihood_ratio + log_prior_ratio;

	double ln_accept_ratio	= log_posterior_ratio + getLnHastingsRatio();

    double lnu = std::log(rng->Uniform(FILE_AND_LINE));
/*	std::cerr << " log_likelihood_ratio = " << log_likelihood_ratio << '\n';
	std::cerr << " log_posterior_ratio = " << log_posterior_ratio << '\n';
	std::cerr << " ln_accept_ratio = " << ln_accept_ratio << '\n'; 
	*/
	if (ln_accept_ratio >= 0.0 || lnu <= ln_accept_ratio)
		{
		p->setLastLnPrior(p->getLastLnPrior() + log_prior_ratio);
		p->setLastLnLike(p->getLastLnLike() + log_likelihood_ratio);
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

