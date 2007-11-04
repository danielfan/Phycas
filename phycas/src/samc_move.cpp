/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
|																			  |
|  This program is free software; you can redistribute it and/or modify		  |
|  it under the terms of the GNU General Public License as published by		  |
|  the Free Software Foundation; either version 2 of the License, or		  |
|  (at your option) any later version.										  |
|																			  |
|  This program is distributed in the hope that it will be useful,			  |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of			  |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the			  |
|  GNU General Public License for more details.								  |
|																			  |
|  You should have received a copy of the GNU General Public License along	  |
|  with this program; if not, write to the Free Software Foundation, Inc.,	  |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.				  |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//#include "phycas/force_include.h"
#include "phycas/src/cipres/CipresDataMatrixHelper.h"
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/samc_move.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/tree_manip.hpp"

//#define KEEP_BRANCHES_LEGAL
//#define MAX_LEGAL_BRLEN 285

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor sets `num_taxa' and `num_nodes_in_fully_resolved_tree', computes the polytomy distribution, and 
|	finally calls SamcMove::reset to re-initialize all other variables.
*/
SamcMove::SamcMove(
  ProbDistShPtr d)	/**< is the distribution to use for generating new edge lengths */
  :MCMCUpdater(),
  term_edge_dist(d),
  view_proposed_move(false)
	{
	num_taxa = 0;
	reset();
#if POLPY_NEWWAY    //SAMC
    goofed = false;
    samc_debug_mode = false;
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void SamcMove::calcPk(
  unsigned leaf_k)	/**< is the leaf we are adding (for extrapolation) or subtracting (for projection) */
	{
	//@POL assumes temperature is infinity, need to rewrite for more reasonable behavior
	pvect.clear();
	unsigned n = tree->GetNNodes() - 1;
	double p = 1.0/(double)n;
	pvect.resize(n, p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
double SamcMove::getPkl(
  unsigned leaf_k,	  /**< is the leaf we are adding (for extrapolation) or subtracting (for projection) */
  TreeNode * nd_l)	  /**< is the leaf we are adding (for extrapolation) or subtracting (for projection) */
	{
	calcPk(leaf_k);
	return pvect[0];	//@POL this will not stand!
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
TreeNode * SamcMove::chooseRandomAttachmentNode(
  unsigned leaf_k)	  /**< is the leaf we are adding (for extrapolation) or subtracting (for projection) */
	{
	calcPk(leaf_k);
	unsigned i = rng->MultinomialDraw(&pvect[0], pvect.size());
	unsigned curr = 0;
	preorder_iterator nd = tree->begin();
	for (++nd; nd != tree->end(); ++nd, ++curr)
		{
		if (curr == i)
			return &(*nd);
		}
	PHYCAS_ASSERT(0);
	return NULL; // avoid non warning
	}

/*----------------------------------------------------------------------------------------------------------------------
|	.
*/
bool SamcMove::extrapolate(
	unsigned leaf_num, /**< */
	double theta_diff,		  /**< the differences in the thetas (this appears in the acceptance probability calculation and is a function of the weights in the SAMC algorithm) */
	double ln_proposal_ratio) /**< the ratio of proposing a projection rather than an extrapolation. */
	{
	const unsigned ninternals_alloced = tree->GetNInternalsAllocated();
	last_move_projection = false;
    if (Tree::gDebugOutput)
    	{
    	std::cerr << "*** extrapolate: before doing anything: " << tree->DebugWalkTree(true, 1) << std::endl; //temporary
	    std::cerr << "*** ninternals_alloced =  " << ninternals_alloced << std::endl; //temporary
		tree->DebugCheckTree(false, true, 2);
		std::cerr << "*** extrapolate: in tree checked" << std::endl; //temporary
		}

    // The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
	// A proposed SamcMove cannot be accepted without changing edge lengths, so it is best to bail out now.
	if (is_fixed)
		return false;

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
	double prev_ln_like = p->getLastLnLike();

    //POL temporary!
    if (samc_debug_mode)
    	{
        likelihood->startTreeViewer(tree, str(boost::format("SamcMove extrapolate (start): leaf_num = %d, prev_ln_like = %f") % leaf_num % prev_ln_like));
        }

    leaf = tree->PopLeafNode();
    PHYCAS_ASSERT(leaf->GetNodeNumber() == leaf_num);
    leaf_sib = chooseRandomAttachmentNode(leaf_num);
    
    if (Tree::gDebugOutput)
        {
        std::cerr << "*** extrapolate:  nodes chosen(but unmodified):" << std::endl; //temporary
        std::cerr << "    leaf: "<< leaf->oneLineDebugReport() << std::endl; //temporary
        std::cerr << "    leaf_sib: "<< leaf_sib->oneLineDebugReport() << std::endl; //temporary
        }
    tree_manipulator.InsertSubtreeIntoEdge(leaf, leaf_sib);
    if (Tree::gDebugOutput)
        {
        std::cerr << "*** extrapolate:  tree structure changed" << std::endl; //temporary
        tree->DebugCheckTree(false, true, 2);
        }
    PHYCAS_ASSERT(ninternals_alloced == tree->GetNInternalsAllocated());
    
    leaf_sib_orig_edgelen = leaf_sib->GetEdgeLen();
    double u = rng->Uniform();
    double new_leaf_sib_edgelen = u*leaf_sib_orig_edgelen;
    double parent_edgelen = leaf_sib_orig_edgelen - new_leaf_sib_edgelen;
    leaf_sib->SetEdgeLen(new_leaf_sib_edgelen);
    new_leaf_sib_parent = parent = leaf_sib->GetParent();
    parent->SetEdgeLen(parent_edgelen);
    double leaf_edgelen = term_edge_dist->Sample();
    leaf->SetEdgeLen(leaf_edgelen);
    
    //POL temporary!
    if (samc_debug_mode)
        {
        leaf_sib->SelectNode();
        parent->SelectNode();
        likelihood->startTreeViewer(tree, str(boost::format("SamcMove extrapolate (after rearrangements, before invalidation), leaf = %d, leaf_sib = %d, parent = %d") % leaf->GetNodeNumber() % leaf_sib->GetNodeNumber() % parent->GetNodeNumber() ));
        parent->UnselectNode();
        leaf_sib->UnselectNode();
        }

    likelihood->useAsLikelihoodRoot(parent);
    likelihood->invalidateAwayFromNode(*parent);
    likelihood->invalidateBothEnds(leaf_sib);
    likelihood->invalidateBothEnds(parent);
    tree->InvalidateNodeCounts();
    PHYCAS_ASSERT(ninternals_alloced == tree->GetNInternalsAllocated());

    //POL temporary!
    if (samc_debug_mode)
        {
        likelihood->startTreeViewer(tree, str(boost::format("SamcMove extrapolate (before calcLnL): leaf = %d, leaf_sib = %d") % leaf->GetNodeNumber() % leaf_sib->GetNodeNumber()));
        }
    
    double curr_ln_like = likelihood->calcLnL(tree);

    //POL temporary!
    if (samc_debug_mode)
        {
        likelihood->startTreeViewer(tree, str(boost::format("SamcMove extrapolate (after calcLnL): curr_ln_like = %f") % curr_ln_like));
        }
    
    double curr_ln_prior = p->calcInternalEdgeLenPriorUnnorm(parent_edgelen);
    curr_ln_prior += p->calcExternalEdgeLenPriorUnnorm(leaf_edgelen);
    double prev_ln_prior = 0.0;
    if (leaf_sib->IsTip())
        {
        prev_ln_prior += p->calcExternalEdgeLenPriorUnnorm(leaf_sib_orig_edgelen);
        curr_ln_prior += p->calcExternalEdgeLenPriorUnnorm(new_leaf_sib_edgelen);
        }
    else
        {
        prev_ln_prior += p->calcInternalEdgeLenPriorUnnorm(leaf_sib_orig_edgelen);
        curr_ln_prior += p->calcInternalEdgeLenPriorUnnorm(new_leaf_sib_edgelen);
        }
    const double prev_posterior = prev_ln_like + prev_ln_prior;
    const double curr_posterior = curr_ln_like + curr_ln_prior;
    double log_pkl_blk_ratio = log(getPkl(leaf_num, leaf_sib)) - log(leaf_sib_orig_edgelen);
    double ln_accept_ratio = theta_diff + curr_posterior - prev_posterior - log_pkl_blk_ratio + ln_proposal_ratio;
    const bool accepted = (ln_accept_ratio >= 0.0 || std::log(rng->Uniform()) <= ln_accept_ratio);
    if (Tree::gDebugOutput)
        {
        std::cerr << "*** extrapolate before revert: " << tree->DebugWalkTree(true, 1) << std::endl; //temporary
        tree->DebugCheckTree(false, true, 2);
        }
    if (accepted)
        {
        p->setLastLnPrior(curr_ln_prior);
        p->setLastLnLike(curr_ln_like);
        accept();
        }
    else
        {
        //curr_ln_like	= p->getLastLnLike();
        //curr_ln_prior	= p->getLastLnPrior();
        revert();
        }
    if (save_debug_info)
        {
        debug_info = str(boost::format("SAMC Extrapolation: deleting %d, curr_post = %f, prev_post = %f, %s") % leaf_num % curr_posterior % prev_posterior % (accepted ? "accepted" : "rejected"));
        }
    PHYCAS_ASSERT(ninternals_alloced == tree->GetNInternalsAllocated());
    return accepted;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   The project move proposes deleting the tip node having node number `leaf_num' from the current tree.
*/
bool SamcMove::project(
  unsigned leaf_num,
  double theta_diff, /* the differences in the thetas (this appears in the acceptance probability calculation and is a function of the weights in the SAMC algorithm) */
  double ln_proposal_ratio) /* the ratio of proposing a projection rather than an extrapolation. */
    {
    last_move_projection = true;
    // The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
    // A proposed SamcMove cannot be accepted without changing edge lengths, so it is best to bail out now.
    if (is_fixed)
        return false;
    
    if (num_taxa == 0)
        {
        throw XLikelihood("Must call finalize() before calling update() for SamcMove");
        }
    
    ChainManagerShPtr p = chain_mgr.lock();
    PHYCAS_ASSERT(p);
    double prev_ln_like = p->getLastLnLike();
    
    leaf = tree->FindTipNode(leaf_num);
    PHYCAS_ASSERT(leaf != NULL);
    orig_edgelen = leaf->GetEdgeLen();
    
    parent = leaf->GetParent();
    PHYCAS_ASSERT(parent != NULL);
    par_orig_edgelen = parent->GetEdgeLen();
    
    PHYCAS_ASSERT(parent->CountChildren() == 2); // generalize when/if we use this with the BushMove
    leaf_sib = leaf->FindNextSib();
    leaf_sib_orig_edgelen = leaf_sib->GetEdgeLen();
    
    new_leaf_sib_parent = leaf_sib->GetParent();
    
    // This DeleteLeaf call will also delete parent and set leaf_sib's edge len
    // 	to the sum of it's edge len and the parent's.
    tree_manipulator.DeleteLeaf(leaf);
    TreeNode *lr = 	leaf_sib->GetParent();
    if (!lr->IsInternal())
        lr = leaf_sib;
    likelihood->useAsLikelihoodRoot(lr);
    
    tree->InvalidateNodeCounts();
    double curr_ln_like = likelihood->calcLnL(tree);
    
    double curr_ln_prior = 0.0;
    double prev_ln_prior = p->calcExternalEdgeLenPriorUnnorm(orig_edgelen);
    prev_ln_prior += p->calcInternalEdgeLenPriorUnnorm(par_orig_edgelen);
    double new_leaf_sib_edgelen = leaf_sib->GetEdgeLen();
    if (leaf_sib->IsTip())
        {
        prev_ln_prior += p->calcExternalEdgeLenPriorUnnorm(leaf_sib_orig_edgelen);
        curr_ln_prior += p->calcExternalEdgeLenPriorUnnorm(new_leaf_sib_edgelen);
        }
    else
        {
        prev_ln_prior += p->calcInternalEdgeLenPriorUnnorm(leaf_sib_orig_edgelen);
        curr_ln_prior += p->calcInternalEdgeLenPriorUnnorm(new_leaf_sib_edgelen);
        }
    double prev_posterior = prev_ln_like + prev_ln_prior;
    double curr_posterior = curr_ln_like + curr_ln_prior;
    double log_pkl_blk_ratio = log(getPkl(leaf_num, leaf_sib)) - log(new_leaf_sib_edgelen);
    double ln_accept_ratio = theta_diff + curr_posterior - prev_posterior + log_pkl_blk_ratio + ln_proposal_ratio;
    const bool accepted = (ln_accept_ratio >= 0.0 || std::log(rng->Uniform()) <= ln_accept_ratio);
    if (accepted)
        {
        p->setLastLnPrior(curr_ln_prior);
        p->setLastLnLike(curr_ln_like);
        accept();
        }
    else
        {
        //curr_ln_like	= p->getLastLnLike();
        //curr_ln_prior	= p->getLastLnPrior();
        revert();
        }
    if (save_debug_info)
        {
        debug_info = str(boost::format("SAMC Projection: deleting %d, curr_post = %f, prev_post = %f, %s") % leaf_num % curr_posterior % prev_posterior % (accepted ? "accepted" : "rejected"));
        }
    return accepted;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or
|	revert(), whichever is appropriate.
*/
bool SamcMove::update()
    {
    return false;
    }

/*--------------------------------------------------------------------------------------------------------------------------
|	Sets `num_taxa' to the number of tips in `tree', computes the polytomy distribution, and sets the number of taxa for the
|	`topo_prior_calculator' object to `num_taxa'.
*/
void SamcMove::finalize()
    {
    num_taxa = tree->GetNTips();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Causes tree to be returned to its state just prior to proposing the move.
*/
void SamcMove::revert()
    {
    if (last_move_projection)
        {
        TreeNode * leafnd = tree->PopLeafNode();
        PHYCAS_ASSERT(leafnd == leaf);
        tree_manipulator.InsertSubtreeIntoEdge(leafnd, leaf_sib);
        leaf_sib->SetEdgeLen(leaf_sib->GetEdgeLen() - par_orig_edgelen);
        leaf->SetEdgeLen(orig_edgelen);
        leaf_sib->GetParent()->SetEdgeLen(par_orig_edgelen);
        likelihood->useAsLikelihoodRoot(new_leaf_sib_parent);
        //likelihood->restoreFromCacheAwayFromNode(*orig_lchild);
        //likelihood->restoreFromCacheParentalOnly(orig_lchild);
        //if (view_proposed_move)
        //	{
        //	orig_par->UnselectNode();
        //	likelihood->startTreeViewer(tree, "Delete edge move REVERTED");
        //	}
        }
    else
        {
        //POL temporary!
        if (samc_debug_mode)
            {
            likelihood->startTreeViewer(tree, str(boost::format("SamcMove extrapolate (before revert): curr_ln_like = %f") % curr_ln_like));
            }
        tree_manipulator.DeleteLeaf(leaf);
        leaf_sib->SetEdgeLen(leaf_sib_orig_edgelen);
        TreeNode * lr = (leaf_sib->IsInternal() ? leaf_sib : leaf_sib->GetParent());
        PHYCAS_ASSERT(lr);
        likelihood->useAsLikelihoodRoot(lr);
        //likelihood->restoreFromCacheAwayFromNode(*lr);
        likelihood->invalidateAwayFromNode(*leaf_sib);
        likelihood->invalidateBothEnds(leaf_sib);
        //POL temporary!
        if (samc_debug_mode)
            {
            likelihood->startTreeViewer(tree, str(boost::format("SamcMove extrapolate (after revert): curr_ln_like = %f") % curr_ln_like));
            }
        }
    tree->InvalidateNodeCounts();
    reset();
    if (Tree::gDebugOutput)
        {
        std::cerr << "*** after revert: " << tree->DebugWalkTree(true, 1) << std::endl; //temporary
        tree->DebugCheckTree(false, true, 2);
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	If `yes' is true, subsequent calls to SamcMove::update will pop up a graphical tree viewer to show the edges
|	affected by the move.
*/
void SamcMove::viewProposedMove(bool yes)
    {
    view_proposed_move = yes;
    }

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void SamcMove::accept()
    {
    if (last_move_projection)
        {
        TreeNode *lr = 	(leaf_sib->IsInternal() ? leaf_sib : leaf_sib->GetParent());
        likelihood->useAsLikelihoodRoot(lr);
        likelihood->discardCacheAwayFromNode(*leaf_sib);
        if (view_proposed_move)
            likelihood->startTreeViewer(tree, "Delete edge move ACCEPTED");
        }
    else
        {
        // Keeping added edge, where orig_par was original polytomous node and orig_lchild was the added node
        likelihood->useAsLikelihoodRoot(leaf_sib->GetParent());
        likelihood->discardCacheAwayFromNode(*leaf_sib);
        likelihood->discardCacheBothEnds(leaf_sib);
        if (view_proposed_move)
            likelihood->startTreeViewer(tree, "extrapolate move ACCEPTED");
        }
    likelihood->invalidateAwayFromNode(*new_leaf_sib_parent);
    reset();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Re-initializes all data members that are recomputed each time ProposeNewState is called.
*/
void SamcMove::reset()
    {
    // Workspace used for computing edge length prior
    // Should always be length one.
    //if (one_edgelen.empty())
    //	  one_edgelen.push_back(0.0);
    leaf_sib = NULL;
    leaf = NULL;
    parent = NULL;
    pvect.clear();
    new_leaf_sib_parent = NULL;
    }
