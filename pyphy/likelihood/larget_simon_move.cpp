#include "phycas/force_include.h"
#include "CipresCommlib/CipresDataMatrixHelper.h"
#include "phycas/rand/probability_distribution.hpp"
#include "pyphy/likelihood/likelihood_models.hpp"
#include "pyphy/likelihood/tree_likelihood.hpp"
#include "pyphy/likelihood/xlikelihood.hpp"
#include "pyphy/likelihood/mcmc_chain_manager.hpp"
#include "pyphy/likelihood/larget_simon_move.hpp"
#include "pyphy/phylogeny/basic_tree.hpp"

#include "boost/format.hpp"

#define KEEP_BRANCHES_LEGAL
#define MAX_LEGAL_BRLEN 285

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
void LargetSimonMove::update()
	{
	// The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
	// A proposed LargetSimonMove cannot be accepted without changing edge lengths, so it is best to bail out now.
	if (is_fixed)
		return;

	ChainManagerShPtr p = chain_mgr.lock();
	assert(p);
	double prev_ln_like = p->getLastLnLike();

	proposeNewState();

	curr_ln_like = likelihood->calcLnL(tree);

	double prev_ln_prior = 0.0;
	if (star_tree_proposal)
		{
		one_edgelen[0]		= orig_edge_len;
		prev_ln_prior		= p->partialEdgeLenPrior(one_edgelen);

		one_edgelen[0]		= orig_node->GetEdgeLen();
		curr_ln_prior		= p->partialEdgeLenPrior(one_edgelen);
		}
	else
		{
		three_edgelens[0] = origX;
		three_edgelens[1] = origY;
		three_edgelens[2] = origZ;
		prev_ln_prior = p->partialEdgeLenPrior(three_edgelens);

		three_edgelens[0] = ndX->GetEdgeLen();
		three_edgelens[1] = ndY->GetEdgeLen();
		three_edgelens[2] = ndZ->GetEdgeLen();
		curr_ln_prior = p->partialEdgeLenPrior(three_edgelens);
		}

	double prev_posterior = prev_ln_like + prev_ln_prior;
	double curr_posterior = curr_ln_like + curr_ln_prior;

	double ln_accept_ratio = curr_posterior - prev_posterior + getLnHastingsRatio() + getLnJacobian();

	if (ln_accept_ratio >= 0.0 || std::log(rng->Uniform(FILE_AND_LINE)) <= ln_accept_ratio)
		{
		p->setLastLnPrior(curr_ln_prior);
		p->setLastLnLike(curr_ln_like);
		accept();
		}
	else
		{
		curr_ln_like	= p->getLastLnLike();
		curr_ln_prior	= p->getLastLnPrior();
		revert();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If the tree has only one internal node (i.e. it is the star tree), then `star_tree_proosal' is set to true and the
|	starTreeProposeNewState function is called. If the tree is not the star tree, then `star_tree_proosal' is set to 
|	false and the defaultProposeNewState function is called.
*/
void LargetSimonMove::proposeNewState()
	{
	if (tree->GetNInternals() == 1)
		{
		starTreeProposeNewState();
		star_tree_proposal = true;
		}
	else
		{
		defaultProposeNewState();
		star_tree_proposal = false;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Chooses a random edge and changes its current length m to a new length m* using the following formula, where `lambda' is
|	a tuning parameter.
|>
|	m* = m*exp(lambda*(r.Uniform(FILE_AND_LINE) - 0.5))
|>
*/
void LargetSimonMove::starTreeProposeNewState()
	{
	// Choose edge randomly.
	//
	unsigned numEdges = tree->GetNNodes() - 1;
	unsigned k = rng->SampleUInt(numEdges);
	unsigned i = 0;
	//@POL this loop is crying out for the for_each algorithm
	for (orig_node = tree->GetFirstPreorder(); orig_node != NULL; orig_node = orig_node->GetNextPreorder())
		{
		// All nodes have an edge associated with them except for the root
		//
		if (!orig_node->IsRoot())
			{
			if (i == k)
				{
				orig_edge_len = orig_node->GetEdgeLen();
				break;
				}
			++i;
			}
		}

	// Modify the edge
	//
	double m		= orig_node->GetEdgeLen();
	double mstar	= m*std::exp(lambda*(rng->Uniform(FILE_AND_LINE) - 0.5));
	orig_node->SetEdgeLen(mstar);

	// Invalidate CLAs to ensure next likelihood calculation will be correct
	orig_node->SelectNode();
	likelihood->useAsLikelihoodRoot(orig_node);
	likelihood->invalidateAwayFromNode(*orig_node);
	likelihood->invalidateBothEnds(orig_node);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a local perturbation using a generalization of the algorithm "LOCAL Without a Molecular Clock" described by 
|	Larget and Simon (1999. Mol. Biol. Evol. 16(6): 750-759). This version allows polytomies (except the case of the 
|	star tree, for which the required three-contiguous-edge segment cannot be identified).
|	
|	     X  c  d
|	      \ | /
|          \|/
|	  a  b  Y
|	   \ | /
|	    \|/
|	     u
|	     |
|	     Z <-- may or may not be the tip at which tree is rooted
|	
|	Pick a random interior node Y whose parent is not the root (i.e. avoid the internal node directly connected to the 
|	tip node at which the tree is rooted. Let u be the parent of y. Let Z be a randomly-chosen child of u (note that in 
|	this case u's parent is considered a "child" of u). In the figure above, we (by chance) chose the parent of u to be 
|	Z, but we could have chosen any of u's real children (except Y). a and b are the other "children" of u, not 
|	including Y. Let X be a randomly chosen child of Y (and here a "child" is really a child). c and d are the other 
|	children of Y.
|	
|	      a   b   c   d
|	       \ /     \ /
|	Z ===== u ===== Y ===== X
|	
|	The path represented by the double line above is either contracted or expanded by a factor m*, where
|	
|	m* = m*exp(lambda*(r.Uniform(FILE_AND_LINE) - 0.5))
|	
|	Then, one of {u, Y} is chosen at random to move. Let's say for illustration that u was chosen. u is moved (along 
|	with a and c) to a random point along the main path from Z to X. If this makes the node u cross over node Y, then 
|	the equivalent of an NNI rearrangement is effected:
|	
|	    X  c  d          X  a  b
|	     \ | /            \ | /    In this case, invalidate CLAs away from u and 
|         \|/              \|/     make u the likelihood root
|	 a  b  Y          c  d  u
|	  \ | /   -->      \ | /
|	   \|/              \|/
|	    u                Y
|	    |                |
|	    Z                Z
|
|	If there is no NNI rearrangement, the move will only require adjusting edge lengths. In this case, invalidate CLAs
|	away from Y and make Y the likelihood root.
*/
void LargetSimonMove::defaultProposeNewState()
	{
	double x, y, z, xstar;

	// Make sure all the necessary shared pointers have been set to something meaningful
	//
	assert(rng);
	assert(tree);
	assert(model);
	assert(likelihood);

	// Begin by resetting all the data members involved with reverting a move
	//
	reset();
	topol_changed = false;

	// Must avoid the "subroot" node (only child of the tip serving as the root), so the number of 
	// acceptable nodes is one fewer than the number of internal nodes
	//
	unsigned numAcceptableNodes = tree->GetNInternals() - 1;

	// Select an internal node whose parent is not the root node to serve as ndY, whose branch will form the middle
	// segment in the path of three contiguous segments to be modified.
	//
	unsigned ypos = rng->SampleUInt(numAcceptableNodes);
	unsigned i = 0;
	for (ndY = tree->GetFirstPreorder(); ndY != NULL; ndY = ndY->GetNextPreorder())
		{
		if (ndY->IsInternal() && !ndY->GetParentConst()->IsRoot())
			{
			if (i == ypos)
				break;
			++i;
			}
		}
	assert(ndY->GetLeftChild() != NULL);
	assert(ndY->GetParentConst() != NULL);
	assert(!ndY->GetParent()->IsRoot());

	// Save ndY's edge length in case revert is needed.
	//
	origY = ndY->GetEdgeLen();

	// Set ndX equal to a randomly-chosen child of ndY
	//
	unsigned ychildren = ndY->CountChildren();
	unsigned which_child = rng->SampleUInt(ychildren);
	unsigned k = 0;
	for (ndX = ndY->GetLeftChild(); ndX != NULL; ndX = ndX->GetRightSib())
		{
		if (k == which_child)
			break;
		++k;
		}
	assert(ndX != NULL);
	origX = ndX->GetEdgeLen();

	// Set ndZ equal to a randomly-chosen child of U = ndY->par, ignoring ndY but treating U->par as if it were a child
	// If U->par is chosen to be ndZ, let ndZ equal U instead because the three nodes (ndX, ndY and ndZ) should be the nodes
	// that manage the relevant edge lengths in the move.
	//
	TreeNode * U = ndY->GetParent();
	assert(U != NULL);
	unsigned uchildren = U->CountChildren();
	which_child = rng->SampleUInt(uchildren);
	if (which_child == 0)
		{
		// Selected "child" is actually U's parent
		//
		ndZ = U;
		origZ = U->GetEdgeLen();

		//	Set ndBase to the deepest affected node
		//
		ndBase = U->GetParent();
		}
	else
		{
		// Selected child is one of U's actual children (but cannot be equal to ndY)
		//
		k = 1;
		for (ndZ = U->GetLeftChild(); ndZ != NULL; ndZ = ndZ->GetRightSib())
			{
			if (ndZ == ndY)
				continue;
			else
				{
				if (k == which_child)
					break;
				++k;
				}
			}
		assert(ndZ != NULL);
		origZ = ndZ->GetEdgeLen();

		//	Set ndBase to the deepest affected node
		//
		ndBase = U;
		}

	m = origX + origY + origZ;
	mstar = m*exp(lambda*(rng->Uniform(FILE_AND_LINE) - 0.5));
	x = origX*mstar/m;
	y = origY*mstar/m;
	z = origZ*mstar/m;

	xstar = rng->Uniform(FILE_AND_LINE)*mstar;

	// Decide whether to move ndY down or node U up
	//
	bool moving_Y = true;
	if (rng->Uniform(FILE_AND_LINE) < 0.5)
		moving_Y = false;
	bool moving_U = !moving_Y;

	if (moving_Y && (xstar <= x + y))
		{
		// No change in topology
		//
		ndX->SetEdgeLen(xstar);
		ndY->SetEdgeLen(x + y - xstar);
		ndZ->SetEdgeLen(z);
		}
	else if (moving_U && (xstar <= y + z))
		{
		// No change in topology
		//
		ndX->SetEdgeLen(x);
		ndY->SetEdgeLen(y + z - xstar);
		ndZ->SetEdgeLen(xstar);
		}
	else
		{
		// Topology changed: figure out which nodes to swap
		//
		if (ndZ != U)
			{
			//cerr << "\n=====> simple <=====\n" << endl;

			// Things are simple when ndZ is U's child and ndX is ndY's child
			//
			swap1 = ndX;
			swap2 = ndZ;

			tree_manipulator.NNISwap(swap1, swap2);
			}
		else
			{
			//cerr << "\n=====> complicated <=====\n" << endl;

			// Things are more complicated when the node to swap is U's parent AND there is
			// the possibility of a polytomy at either U or ndY. If we could assume that polytomies
			// were impossible, we could make our lives easier by just swapping U's other "child"
			// (i.e. the one that is not U's parent and also not ndY) and ndY's other child (i.e. the one that
			// is not ndX); however, if there are possibly more children of either node U or ndY, then
			// this will not work.
			//
			// Remember that ndZ is actually the same node as U in this case. If we later need to revert
			// this move, we know to avoid using NNISwap if ndZ->GetParent() == swap2. This is the reason
			// that we go to the trouble of setting swap2 here, even though it is not needed by NNISwapSpecial.
			//
			swap1 = ndX;
			swap2 = ndZ->GetParent();

			tree_manipulator.NNISwapSpecial(swap1);
			}

		if (moving_Y)
			{
			ndX->SetEdgeLen(x + y);
			ndY->SetEdgeLen(xstar - x - y);
			ndZ->SetEdgeLen(x + y + z - xstar);
			}
		else
			{
			ndX->SetEdgeLen(x + y + z - xstar);
			ndY->SetEdgeLen(xstar - y - z);
			ndZ->SetEdgeLen(y + z);
			}

		assert(ndX->GetEdgeLen() > 0.0);
		assert(ndY->GetEdgeLen() > 0.0);
		assert(ndZ->GetEdgeLen() > 0.0);

		topol_changed = true;
		}

	ndX->SelectNode();
	ndY->SelectNode();
	ndZ->SelectNode();
	likelihood->useAsLikelihoodRoot(ndY);
	likelihood->invalidateAwayFromNode(*ndY);
	likelihood->invalidateBothEnds(ndY);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reverses move made in proposeNewState. Assumes ndX, ndY, and ndZ are non-NULL, which will be true if proposeNewState
|	was just called.
*/
void LargetSimonMove::revert()
	{
	if (star_tree_proposal)
		{
		orig_node->SetEdgeLen(orig_edge_len);

		likelihood->useAsLikelihoodRoot(orig_node);
		likelihood->restoreFromCacheAwayFromNode(*orig_node);
		likelihood->invalidateBothEndsDiscardCache(orig_node);

		orig_node->UnselectNode();
		}
	else
		{
		assert(ndX != NULL);
		assert(ndY != NULL);
		assert(ndZ != NULL);
		assert(topol_changed ? (swap1 != NULL && swap2 != NULL) : (swap1 == NULL && swap2 == NULL));

		if (topol_changed)
			{
			if (swap2 == ndZ)
				{
				// If swap2 equals ndZ, then swap2 was a child of ndBase and we were able to use
				// the standard NNISwap function to swap the two nodes
				//
				tree_manipulator.NNISwap(swap1, swap2);
				}
			else 
				{
				// If swap2 is ndZ's parent, then swap2 is ndBase (i.e. it is the "child" node below the
				// lower of the two adjacent internal nodes involved in the swap) and we had to use the
				// NNISwapSpecial function to perform the rearrangment
				//
				tree_manipulator.NNISwapSpecial(swap1);
				}
			}

		ndX->SetEdgeLen(origX);
		ndY->SetEdgeLen(origY);
		ndZ->SetEdgeLen(origZ);

		likelihood->useAsLikelihoodRoot(ndY);
		likelihood->restoreFromCacheAwayFromNode(*ndY);
		likelihood->invalidateBothEndsDiscardCache(ndY);

		ndX->UnselectNode();
		ndY->UnselectNode();
		ndZ->UnselectNode();
		}

	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void LargetSimonMove::accept()
	{
	likelihood->useAsLikelihoodRoot(ndY);
	likelihood->discardCacheAwayFromNode(*ndY);
	likelihood->discardCacheBothEnds(ndY);

	ndX->UnselectNode();
	ndY->UnselectNode();
	ndZ->UnselectNode();

	reset();
	}

}	// namespace phycas
