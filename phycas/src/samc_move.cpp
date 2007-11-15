/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2007 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
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

#include "phycas/src/cipres/CipresDataMatrixHelper.h"
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/samc_move.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/tree_manip.hpp"

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor sets `num_taxa' and `num_nodes_in_fully_resolved_tree', computes the polytomy distribution, and 
|	finally calls SamcMove::reset to re-initialize all other variables.
*/
SamcMove::SamcMove(
  ProbDistShPtr d)	/**< is the distribution to use for generating new edge lengths */
  :MCMCUpdater(),
  term_edge_dist(d),
  temperature_threshhold(100.0)
	{
	num_taxa = 0;
	reset();
#if POLPY_NEWWAY
    prev_ln_like = 0.0;
    infinite_temperature = false;
    temperature = 0.6;
#endif
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|   Sets the maximum number of taxa. This is used to ensure that Split objects in the tree have enough bits to 
|   accommodate any combination of taxa that could appear in the tree at the same time.
*/
void SamcMove::setNTax(
  unsigned ntax)   /**< is the maximum number of taxa that could appear in the tree */
    {
    num_taxa = ntax;
    }
#endif

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|   Returns the current value of the data member `num_taxa'.
*/
unsigned SamcMove::getNTax()
    {
    return num_taxa;
    }
#endif

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Computes the least squares length (v_i) of an edge connecting `leaf_k' to the edge subtending node `nd'. The formula
|   used is eq. 4 from Rzhetsky and Nei (1993. MBE 10:1074). This is used by SamcMove::calcPk to compute the probability
|   of selecting a node in the tree (or rather an edge in the tree) to which the new leaf will be attached in an 
|   extrapolate proposal.
*/
double SamcMove::calcLeafEdgeLen(
  unsigned leaf_k,  /**< is the leaf we are adding in an extrapolation move */
  TreeNode * nd)	/**< is the putative attachment node */
	{
    // Start by getting the set of all nodes above nd in the tree
    std::vector<unsigned> nodes_above;
    nd->GetSplit().GetOnListImpl(nodes_above);
    PHYCAS_ASSERT(nodes_above.size() > 0);
    double na = (double)nodes_above.size();

    //std::cerr << "  leaf_k = " << leaf_k << "\n";
    //std::cerr << "   " << na << " nodes above: ";
    //std::copy(nodes_above.begin(), nodes_above.end(), std::ostream_iterator<double>(std::cerr," "));
    //std::cerr << std::endl;

    // Now get the set of all nodes below nd in the tree
    std::vector<unsigned> nodes_below;
    nd->GetSplit().GetOffListImpl(nodes_below);
    PHYCAS_ASSERT(nodes_below.size() > 0);
    double nb = (double)nodes_below.size();

    //std::cerr << "   " << nb << " nodes below: ";
    //std::copy(nodes_below.begin(), nodes_below.end(), std::ostream_iterator<double>(std::cerr," "));
    //std::cerr << std::endl;

    // Calculate sum of all nabove distances between leaf_k and each node in nodes_above vector
    // At same time, compute ntotal distances between each node in nodes_above vector with each node in
    // the nodes_below vector
    double d_xa = 0.0;
    double d_ab = 0.0;
    for (std::vector<unsigned>::const_iterator ait = nodes_above.begin(); ait != nodes_above.end(); ++ait)
        {
        d_xa += distance_matrix[leaf_k][*ait];
        for (std::vector<unsigned>::const_iterator bit = nodes_below.begin(); bit != nodes_below.end(); ++bit)
            {
            d_ab += distance_matrix[*ait][*bit];
            }
        }

    // Calculate sum of all nbelow distances between leaf_k and each node in nodes_below vector
    double d_xb = 0.0;
    for (std::vector<unsigned>::const_iterator bit = nodes_below.begin(); bit != nodes_below.end(); ++bit)
        {
        d_xb += distance_matrix[leaf_k][*bit];
        }

    // Now in a position to compute the least squares length of the edge subtending leaf_k if it were
    // added to the tree on the edge below node nd
    double d = (d_xa/na + d_xb/nb - d_ab/(na*nb))/2.0;

    //std::cerr << "   d_xa = " << d_xa << std::endl;
    //std::cerr << "   d_xb = " << d_xb << std::endl;
    //std::cerr << "   d_ab = " << d_ab << std::endl;
    //std::cerr << "   d    = " << d << std::endl;

    return d;
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the length (v_i) of the edge connecting `leak_k' to every existing edge i in the tree. The formula used is
|   eq. 4 from Rzhetsky and Nei (1993. MBE 10:1074). A better method would be to compute the minimum evolution score
|   resulting from attaching `leaf_k' to every possible edge in the existing tree, but would be much more expensive so
|   the approximation is used instead. For each i, this function then computes the relative probability of insertion 
|   p*_i = exp{-v_i/tau_s}, where tau_s is the temperature parameter. The normalized values p_i = p*_i/sum_i(p*_i) are
|   stored in the vector `pvect' for use by SamcMove::chooseRandomAttachmentNode in selecting a node in the tree (or 
|   rather an edge in the tree) to which the new `leaf_k' will be attached in an extrapolate proposal.
*/
void SamcMove::calcPk(
  unsigned leaf_k)	/**< is the leaf we are adding (for extrapolation) or subtracting (for projection) */
	{
	pvect.clear();
#if POLPY_NEWWAY
    tree->RecalcAllSplits(num_taxa);
    if (infinite_temperature)
        {
	    // Efficient version when temperature is infinity
	    unsigned n = tree->GetNNodes() - 1;
	    double p = 1.0/(double)n;
    	pvect.resize(n, p);
        }
    else
        {
        //std::cerr << "\nCURRENT TREE: " << tree->MakeNewick() << std::endl;
	    TreeNode * nd = tree->GetFirstPreorder();
	    for (nd = nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		    {
            //std::cerr << "***** calcLeafEdgeLen: leaf_k = " << leaf_k << ", nd = " << nd->GetNodeNumber() << std::endl;
            double d = calcLeafEdgeLen(leaf_k, nd);
            double p = std::exp(-d/temperature);
            pvect.push_back(p);
		    }

        //std::cerr << "\npvect before normalizing:" << std::endl;
        //std::copy(pvect.begin(), pvect.end(), std::ostream_iterator<double>(std::cerr,"|"));

        PHYCAS_ASSERT(pvect.size() == tree->GetNNodes() - 1);
    	double total = std::accumulate(pvect.begin(), pvect.end(), 0.0, std::plus<double>());
        std::for_each(pvect.begin(), pvect.end(), boost::lambda::_1 /= total);

        //std::cerr << "\n\npvect after normalizing:" << std::endl;
        //std::copy(pvect.begin(), pvect.end(), std::ostream_iterator<double>(std::cerr,"|"));
        //std::cerr << std::endl;
        }
#else
	// Efficient version when temperature is infinity
	unsigned n = tree->GetNNodes() - 1;
	double p = 1.0/(double)n;
	pvect.resize(n, p);
#endif
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
|	Proposes adding next taxon in the ladder. Need to be generalized so that more than one taxon can be added at once.
*/
bool SamcMove::extrapolate(
	unsigned leaf_num,          /**< is the index of the taxon to be added */
	double theta_diff,		    /**< is the difference in the thetas (this appears in the acceptance probability calculation and is a function of the weights in the SAMC algorithm) */
	double ln_proposal_ratio)   /**< is the ratio of proposing a projection rather than an extrapolation, also needed in the acceptance probability. */
	{
	const unsigned ninternals_alloced = tree->GetNInternalsAllocated();
	last_move_projection = false;

    // The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
	// A proposed SamcMove cannot be accepted without changing edge lengths, so it is best to bail out now.
	if (is_fixed)
		return false;

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
	prev_ln_like = p->getLastLnLike();

    leaf = tree->PopLeafNode();
    PHYCAS_ASSERT(leaf->GetNodeNumber() == leaf_num);
    leaf_sib = chooseRandomAttachmentNode(leaf_num);
    
    // Examples of debugging tools:
    // if (tree->debugOutput)
    //   {
    //   std::cerr << "leaf: " << leaf->oneLineDebugReport() << std::endl;
    //   leaf_sib->SelectNode();
    //   likelihood->startTreeViewer(tree, str(boost::format("SamcMove extrapolate (start): leaf_num = %d, prev_ln_like = %f") % leaf_num % prev_ln_like));
    //   leaf_sib->UnselectNode();
    //   }

    tree_manipulator.InsertSubtreeIntoEdge(leaf, leaf_sib);
    PHYCAS_ASSERT(ninternals_alloced == tree->GetNInternalsAllocated());
    
    leaf_sib_orig_edgelen = leaf_sib->GetEdgeLen();
    double u = rng->Uniform(FILE_AND_LINE);
    double new_leaf_sib_edgelen = u*leaf_sib_orig_edgelen;
    double parent_edgelen = leaf_sib_orig_edgelen - new_leaf_sib_edgelen;
    leaf_sib->SetEdgeLen(new_leaf_sib_edgelen);
    new_leaf_sib_parent = parent = leaf_sib->GetParent();
    parent->SetEdgeLen(parent_edgelen);
    double leaf_edgelen = term_edge_dist->Sample();
    leaf->SetEdgeLen(leaf_edgelen);
    
    likelihood->useAsLikelihoodRoot(parent);
    likelihood->invalidateAwayFromNode(*parent);
    likelihood->invalidateBothEnds(leaf_sib);
    likelihood->invalidateBothEndsDiscardCache(parent); // wouldn't be necessary if nodes retrieved from storage were guaranteed to be clean
    tree->InvalidateNodeCounts();
    PHYCAS_ASSERT(ninternals_alloced == tree->GetNInternalsAllocated());

    double curr_ln_like = likelihood->calcLnL(tree);

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
    const bool accepted = (ln_accept_ratio >= 0.0 || std::log(rng->Uniform(FILE_AND_LINE)) <= ln_accept_ratio);
    if (accepted)
        {
        p->setLastLnPrior(curr_ln_prior);
        p->setLastLnLike(curr_ln_like);
        accept();
        }
    else
        {
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
  unsigned leaf_num,        /**< is the number of the node to propose deleting */
  double theta_diff,        /**< the differences in the thetas (this appears in the acceptance probability calculation and is a function of the weights in the SAMC algorithm) */
  double ln_proposal_ratio) /**< the ratio of proposing a projection rather than an extrapolation. */
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
    prev_ln_like = p->getLastLnLike();
    
    // The next line assumes that tip serving as root will never be chosen. This is (currently) guaranteed 
    // by the fact that the tip node serving as the root is one of the original four taxa used to build
    // the starting tree for SAMC and is never removed because the lowest level is the four-taxon level.
    leaf = tree->FindTipNode(leaf_num);
    PHYCAS_ASSERT(leaf != NULL);
    orig_edgelen = leaf->GetEdgeLen();
    
    parent = leaf->GetParent();
    PHYCAS_ASSERT(parent != NULL);
    par_orig_edgelen = parent->GetEdgeLen();
    
    PHYCAS_ASSERT(parent->CountChildren() == 2); // generalize when/if we use this with the BushMove
    leaf_sib = leaf->FindNextSib();
    leaf_sib_orig_edgelen = leaf_sib->GetEdgeLen();
    
    // This DeleteLeaf call will also delete parent and will add parent's edge len to leaf_sib's edge len
    tree_manipulator.DeleteLeaf(leaf);
    new_leaf_sib_parent = leaf_sib->GetParent();
    if (!new_leaf_sib_parent->IsInternal())
        new_leaf_sib_parent = leaf_sib;
    likelihood->useAsLikelihoodRoot(new_leaf_sib_parent);
    likelihood->invalidateAwayFromNode(*new_leaf_sib_parent);
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
    const bool accepted = (ln_accept_ratio >= 0.0 || std::log(rng->Uniform(FILE_AND_LINE)) <= ln_accept_ratio);
    if (accepted)
        {
        p->setLastLnPrior(curr_ln_prior);
        p->setLastLnLike(curr_ln_like);
        accept();
        }
    else
        {
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
//void SamcMove::finalize()
//    {
//    num_taxa = tree->GetNTips();
//    }

/*----------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Causes tree to be returned to its state just prior to proposing the move.
*/
void SamcMove::revert()
    {
    if (last_move_projection)
        {
        // Restore from cache before performing the reverse move because we originally invalidated after the forward move
        likelihood->restoreFromCacheAwayFromNode(*new_leaf_sib_parent);

        // Reverse move
        TreeNode * leafnd = tree->PopLeafNode();
        PHYCAS_ASSERT(leafnd == leaf);
        tree_manipulator.InsertSubtreeIntoEdge(leafnd, leaf_sib);
        leaf_sib->SetEdgeLen(leaf_sib->GetEdgeLen() - par_orig_edgelen);
        leaf->SetEdgeLen(orig_edgelen);
        leaf_sib->GetParent()->SetEdgeLen(par_orig_edgelen);

        // Be sure to invalidate both ends of node used as likelihood root in case one or both of these CLAs were 
        // invalid before the forward move and are now valid as a result of computing the likelihood after the forward move
        likelihood->invalidateBothEndsDiscardCache(new_leaf_sib_parent);
        }
    else
        {
        likelihood->restoreFromCacheAwayFromNode(*parent);
        likelihood->restoreFromCacheBothEnds(leaf_sib);
        tree_manipulator.DeleteLeaf(leaf);
        leaf_sib->SetEdgeLen(leaf_sib_orig_edgelen);
        TreeNode * lr = (leaf_sib->IsInternal() ? leaf_sib : leaf_sib->GetParent());
        PHYCAS_ASSERT(lr);
        likelihood->useAsLikelihoodRoot(lr);
        }
    tree->InvalidateNodeCounts();
    reset();
    }

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void SamcMove::accept()
    {
    if (last_move_projection)
        {
        likelihood->discardCacheBothEnds(new_leaf_sib_parent);
        likelihood->discardCacheAwayFromNode(*new_leaf_sib_parent); // POL was leaf_sib but that was incorrect(?)
        }
    else
        {
        // Keeping added edge, where orig_par was original polytomous node and orig_lchild was the added node
        likelihood->useAsLikelihoodRoot(leaf_sib->GetParent());
        likelihood->discardCacheAwayFromNode(*parent);
        likelihood->discardCacheBothEnds(leaf_sib);
        }
    reset();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Re-initializes all data members that are recomputed each time ProposeNewState is called.
*/
void SamcMove::reset()
    {
    leaf_sib = NULL;
    leaf = NULL;
    parent = NULL;
    pvect.clear();
    new_leaf_sib_parent = NULL;
    }

#if POLPY_NEWWAY    //SAMC
/*----------------------------------------------------------------------------------------------------------------------
|	Saves one row of the distance matrix to the data member `distance_matrix'.
*/
void SamcMove::setDistanceMatrixRow(
  const DistanceMatrixRow & row)    /**< is the row to save */
    {
    DistanceMatrixRow r;
    r.resize(row.size(), 0);
    std::copy(row.begin(), row.end(), r.begin());
    distance_matrix.push_back(r);
    }
#endif

#if POLPY_NEWWAY    //SAMC
/*----------------------------------------------------------------------------------------------------------------------
|	Setter for data member `temperature'.
*/
void SamcMove::setTemperature(
  double temp)  /**< is the new temperature */
    {
    PHYCAS_ASSERT(temp > 0.0);
    temperature = temp;
    if (temperature > temperature_threshhold)
        infinite_temperature = true;
    else
        infinite_temperature = false;
    }
#endif

#if POLPY_NEWWAY    //SAMC
/*----------------------------------------------------------------------------------------------------------------------
|	Accessor for data member `temperature'.
*/
double SamcMove::getTemperature() const
    {
    return temperature;
    }
#endif
