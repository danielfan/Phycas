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
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/tip_data.hpp"
#include "phycas/src/internal_data.hpp"
#include "phycas/src/sim_data.hpp"
#include "phycas/src/phycas_string.hpp"
#include "phycas/src/basic_lot.hpp"
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <numeric>
#include "phycas/src/edge_iterators.hpp"

static int8_t codon_state_codes[] =
	{
	0,  // 0 AAA
	1,  // 1 AAC
	2,  // 2 AAG
	3,  // 3 AAT
	4,  // 4 ACA
	5,  // 5 ACC
	6,  // 6 ACG
	7,  // 7 ACT
	8,  // 8 AGA
	9,  // 9 AGC
	10, // 10 AGG
	11, // 11 AGT
	12, // 12 ATA
	13, // 13 ATC
	14, // 14 ATG
	15, // 15 ATT
	16, // 16 CAA
	17, // 17 CAC
	18, // 18 CAG
	19, // 19 CAT
	20, // 20 CCA
	21, // 21 CCC
	22, // 22 CCG
	23, // 23 CCT
	24, // 24 CGA
	25, // 25 CGC
	26, // 26 CGG
	27, // 27 CGT
	28, // 28 CTA
	29, // 29 CTC
	30, // 30 CTG
	31, // 31 CTT
	32, // 32 GAA
    33, // 33 GAC
	34, // 34 GAG
	35, // 35 GAT
	36, // 36 GCA
	37, // 37 GCC
	38, // 38 GCG
	39, // 39 GCT
	40, // 40 GGA
	41, // 41 GGC
	42, // 42 GGG
	43, // 43 GGT
	44, // 44 GTA
	45, // 45 GTC
	46, // 46 GTG
	47, // 47 GTT
	61, // 48 TAA stop
	48, // 49 TAC
	61, // 50 TAG stop
	49, // 51 TAT
	50, // 52 TCA
	51, // 53 TCC
	52, // 54 TCG
	53, // 55 TCT
	61, // 56 TGA stop
	54, // 57 TGC
	55, // 58 TGG
	56, // 59 TGT
	57, // 60 TTA
	58, // 61 TTC
	59, // 62 TTG
	60  // 63 TTT
	};

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Updates the conditional likelihood array for `nd' in a direction away from `avoid'. All adjacent nodes (other than 
|	`avoid') are assumed to be valid.
*/
void TreeLikelihood::refreshCLA(TreeNode & nd, const TreeNode * avoid)
	{
	if (nd.IsTip())
		return;

	PHYCAS_ASSERT(avoid != NULL);
	//@MTH-NESCENT this const_cast should be safe because we aren't relying on const-ness of CLA's but it needs to be revisited
	//@POL-NESCENT Mark, if we are going to immediately cast away the const on avoid, why do we need to make it const in the first place?
	CondLikelihoodShPtr ndCondLike = getCondLikePtr(&nd, const_cast<TreeNode *>(avoid)); 
	TreeNode * parent = nd.GetParent();
	TreeNode * lChild = nd.GetLeftChild();
	PHYCAS_ASSERT(parent != NULL);
	PHYCAS_ASSERT(lChild != NULL);

	// The first neighbor can either be the parent or the leftmost child. The second neighbor must be a child, but
	// which child depends on the first neighbor. Third and subsequent neighbors are always next sibs.
	TreeNode * firstNeighbor = NULL;
	TreeNode * secondNeighbor = NULL;
	double firstEdgeLen = 0.0;
	const bool movingTowardLeaves = !(parent == avoid);
	if (movingTowardLeaves)
		{
		// A child of nd is closer to the likelihood root than nd
		firstNeighbor = parent;
		firstEdgeLen = nd.GetEdgeLen();
		secondNeighbor = (lChild == avoid ? lChild->GetRightSib() : lChild);
		}
	else
		{
		// The parent of nd is closer to the likelihood root than nd
		firstNeighbor = lChild;
		firstEdgeLen = lChild->GetEdgeLen();
		secondNeighbor = lChild->GetRightSib();
		}
	
	if (firstNeighbor->IsTip())
		{
		TipData & firstTD = *(firstNeighbor->GetTipData());
		calcPMatTranspose(firstTD.getTransposedPMatrices(), firstTD.getConstStateListPos(), firstEdgeLen);
		if (secondNeighbor->IsTip())
			{
			// 1. both neighbors are tips
			TipData & secondTD = *(secondNeighbor->GetTipData());
			calcPMatTranspose(secondTD.getTransposedPMatrices(), secondTD.getConstStateListPos(), secondNeighbor->GetEdgeLen());
			calcCLATwoTips(*ndCondLike, firstTD, secondTD);
			}
		else
			{
			// 2. first neighbor is a tip, but second is an internal node
			InternalData & secondID	= *(secondNeighbor->GetInternalData());
			calcPMat(secondID.getPMatrices(), secondNeighbor->GetEdgeLen());
			CondLikelihoodShPtr secCL = getCondLikePtr(secondNeighbor, &nd);
			ConstPMatrices secPMat = secondID.getConstPMatrices();
			calcCLAOneTip(*ndCondLike, firstTD, secPMat, *secCL);
			}
		}
	else
		{
		InternalData & firstID = *(firstNeighbor->GetInternalData());
		calcPMat(firstID.getPMatrices(), firstEdgeLen);
		const CondLikelihood & firCL = *getCondLikePtr(firstNeighbor, &nd);
		ConstPMatrices firPMat = firstID.getConstPMatrices();	
		if (secondNeighbor->IsTip())
			{
			// 3. first neighbor internal node, but second is a tip
			TipData & secondTD = *(secondNeighbor->GetTipData());
			calcPMatTranspose(secondTD.getTransposedPMatrices(), secondTD.getConstStateListPos(), secondNeighbor->GetEdgeLen());
			calcCLAOneTip(*ndCondLike, secondTD, firPMat, firCL);
			}
		else
			{
			// 4. both neighbors are internal nodes
			InternalData & secondID	= *(secondNeighbor->GetInternalData());
			calcPMat(secondID.getPMatrices(), secondNeighbor->GetEdgeLen());
			const CondLikelihood & secCL = *getCondLikePtr(secondNeighbor, &nd);
			ConstPMatrices secPMat = secondID.getConstPMatrices();
			calcCLANoTips(*ndCondLike, firPMat, firCL, secPMat, secCL);
			}
		}

	// Deal with possible polytomy in which secondNeighbor has siblings
	for (TreeNode * currNd = secondNeighbor->GetRightSib(); currNd != NULL; currNd = currNd->GetRightSib())
		{
		if (currNd != avoid)
			{
			if (currNd->IsTip())
				{
				TipData & currTD = *(currNd->GetTipData());
				calcPMatTranspose(currTD.getTransposedPMatrices(), currTD.getConstStateListPos(), currNd->GetEdgeLen());
				conditionOnAdditionalTip(*ndCondLike, currTD);
				}
			else
				{
				InternalData & currID = *(currNd->GetInternalData());
				calcPMat(currID.getPMatrices(), currNd->GetEdgeLen());
				const CondLikelihood & currCL = *getCondLikePtr(currNd, &nd);
				ConstPMatrices currPMat = currID.getConstPMatrices();
				conditionOnAdditionalInternal(*ndCondLike, currPMat, currCL);
				}
			}
		}

#if 0 && POLPY_NEWWAY
    // Turn this section on for debugging purposes only! Walks through conditional likelihood arrays
    // for this node looking for any conditional likelihood that is negative. Negative cond. likes
    // have shown up in the past (18 Oct 2007) for the GTR model, which results in a NaN (represented
    // as -1.#IND in the VC compiler when the log of the site likelihood is taken.
	LikeFltType * cla = ndCondLike->getCLA();
	for (unsigned i = 0; i < num_rates; ++i)
		{
	    for (unsigned j = 0; j < num_patterns; ++j)
		    {
	        for (unsigned k = 0; k < num_states; ++k)
		        {
                if (*cla++ < 0.0)
                    {
                    std::cerr << "*** Doh! cla negative for rate = " << i << ", pattern = " << j << ", state = " << k << " ***" << std::endl;
                    std::exit(0);
                    }
		        }
	        }
        }

#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the conditional likelihood of refNd is up to date for calculations centered at some effective root 
|	node (neighborCloserToEffectiveRoot will be a node adjacent to refNd, but closer than refNd to the the effective 
|	root). Called in a context in which neighborCloserToEffectiveRoot is requesting that all of its neighbors update 
|	their likelihood temporaries.
*/
bool TreeLikelihood::isValid(const TreeNode * refNd, const TreeNode * neighborCloserToEffectiveRoot)
	{
	PHYCAS_ASSERT(refNd != NULL);
	PHYCAS_ASSERT(neighborCloserToEffectiveRoot != NULL);
	PHYCAS_ASSERT(refNd != neighborCloserToEffectiveRoot);
	if (refNd->IsTip())
		{
		// tip nodes always return true because they must be the child of neighborCloserToEffectiveRoot
		// and they do not hold filial CLAs
		return true;
		}
	else
		{
		// refNd is internal
		if (refNd->GetParentConst() == neighborCloserToEffectiveRoot)
			{
			// refNd is the child of neighborCloserToEffectiveRoot
			// If refNd has a filial CLA, then return true because that means refNd is valid
			const InternalData * id = refNd->GetInternalData();
			return (id->childWorkingCLA);			
			}
		else
			{
			// neighborCloserToEffectiveRoot is the child of refNd
			// If neighborCloserToEffectiveRoot has a parental CLA, then return true because 
			// that means refNd is valid
			const InternalData * id = neighborCloserToEffectiveRoot->GetInternalData();
			return (id->parWorkingCLA);			
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Unconditionally invalidates parental (and, if ref_nd is internal, filial) conditional likelihood arrays and also
|	removes cached conditional likelihood arrays if present. This function always returns false so it can be used in 
|	conjunction with effective_postorder_edge_iterator to invalidate every CLA in the entire tree and ensure that there
|	are also no cached CLAs as well.
*/
bool TreeLikelihood::invalidateBothEndsDiscardCache(TreeNode * ref_nd, TreeNode * /* unused */)
	{
	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();
		if (td != NULL)
			{
			// Invalidate the parental CLAs if they exist
			if (td->parWorkingCLA)
				{
				cla_pool.putCondLikelihood(td->parWorkingCLA);
				td->parWorkingCLA.reset();
				}
			// Remove cached parental CLAs if they exist
			if (td->parCachedCLA)
				{
				cla_pool.putCondLikelihood(td->parCachedCLA);
				td->parCachedCLA.reset();
				}
			}
		}
	else
		{
		// Internal nodes have both parental and filial CLAs
		InternalData * id = ref_nd->GetInternalData();
		if (id != NULL)
			{
			// Invalidate the parental CLAs if they exist
			if (id->parWorkingCLA)
				{
				cla_pool.putCondLikelihood(id->parWorkingCLA);
				id->parWorkingCLA.reset();
				}
			// Remove cached parental CLAs if they exist
			if (id->parCachedCLA)
				{
				cla_pool.putCondLikelihood(id->parCachedCLA);
				id->parCachedCLA.reset();
				}

			// Invalidate the filial CLAs if they exist
			if (id->childWorkingCLA)
				{
				cla_pool.putCondLikelihood(id->childWorkingCLA);
				id->childWorkingCLA.reset();
				}
			// Remove cached filial CLAs if they exist
			if (id->childCachedCLA)
				{
				cla_pool.putCondLikelihood(id->childCachedCLA);
				id->childCachedCLA.reset();
				}
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Unconditionally invalidate parental (and, if ref_nd is internal, filial) conditional likelihood arrays. Always 
|	returns false so it can be used in conjunction with effective_postorder_edge_iterator to invalidate every CLA in
|	the entire tree.
*/
bool TreeLikelihood::invalidateBothEnds(TreeNode * ref_nd, TreeNode * /* unused */)
	{
	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();

		// Invalidate the parental CLAs if they exist
		if (td->parWorkingCLA)
			{
			if (td->parCachedCLA)
				cla_pool.putCondLikelihood(td->parCachedCLA);
			td->parCachedCLA = td->parWorkingCLA;
			td->parWorkingCLA.reset();
			}
		}
	else
		{
		// Internal nodes have both parental and filial CLAs
		InternalData * id = ref_nd->GetInternalData();

		// Invalidate the parental CLAs if they exist
		if (id->parWorkingCLA)
			{
			if (id->parCachedCLA)
				cla_pool.putCondLikelihood(id->parCachedCLA);
			id->parCachedCLA = id->parWorkingCLA;
			id->parWorkingCLA.reset();
			}

		// Invalidate the filial CLAs if they exist
		if (id->childWorkingCLA)
			{
			if (id->childCachedCLA)
				cla_pool.putCondLikelihood(id->childCachedCLA);
			id->childCachedCLA = id->childWorkingCLA;
			id->childWorkingCLA.reset();
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function used in conjunction with effective_postorder_edge_iterator to invalidate the appropriate conditional 
|	likelihood arrays (CLAs) of all nodes starting with a focal node. For nodes that are descendants of the focal node
|	(descendant means that a node can be found using only left child and right sib pointers starting from the focal
|	node), it is the parental CLAs that are invalidated. For a node that is an ancestor of the focal node, it is the 
|	filial CLAs that are invalidated. For all other nodes (e.g. on independent lineages derived from an ancestral node), 
|	it is the parental CLAs that are invalidated. This pattern of invalidation ensures that the likelihood will be 
|	correctly computed using any node in the tree as the likelihood root. For any likelihood root n, it now appears as 
|	though an edge length just on the other side of the focal node from n has been changed, necessitating the 
|	recalculation of all CLAs starting from that point back toward the likelihood root. The return value indicates
|	whether or not the iterator should continue into the subtree defined by `refNd' (away from 
|	`neighborCloserToEffectiveRoot'). This function always returns false because the goal is to invalidate every node
|	that needs to be invalidated.
*/
bool TreeLikelihood::invalidateNode(TreeNode * ref_nd, TreeNode * neighbor_closer_to_likelihood_root)
	{
	if (ref_nd->IsTip() && !ref_nd->IsTipRoot())
		{
		TipData * td = ref_nd->GetTipData();
		if (!td->parWorkingCLA)
			return false;
		if (td->parCachedCLA)
			cla_pool.putCondLikelihood(td->parCachedCLA);
		td->parCachedCLA = td->parWorkingCLA;
		td->parWorkingCLA.reset();
		}
	else
		{
		if (ref_nd->GetParent() == neighbor_closer_to_likelihood_root)
			{
			// ref_nd is the actual child
			InternalData * id = ref_nd->GetInternalData();
			if (!id->parWorkingCLA)
				return false;
			if (id->parCachedCLA)
				cla_pool.putCondLikelihood(id->parCachedCLA);
			id->parCachedCLA = id->parWorkingCLA;
			id->parWorkingCLA.reset();
			}
		else
			{
			// neighbor_closer_to_likelihood_root is the actual child
			if (neighbor_closer_to_likelihood_root->IsTip())
				return false;
			// Either ref_nd is not a tip, or it is the tip serving as the root node
			InternalData * id = neighbor_closer_to_likelihood_root->GetInternalData();
			if (!id->childWorkingCLA)
				return false;
			if (id->childCachedCLA)
				cla_pool.putCondLikelihood(id->childCachedCLA);
			id->childCachedCLA = id->childWorkingCLA;
			id->childWorkingCLA.reset();
			}
		}
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Invalidates all conditional likelihood arrays (CLAs) pointing toward the supplied `focal_node'. That is, the
|	parental CLA of a child of `focal_node' would be invalidated, as would be the filial CLA for the parent of the
|	`focal_node'. The invalidation propagates throughout the tree from `focal_node' so that regardless of the node
|	serving as the likelihood root, the correct CLAs will be recalculated the next time the log-likelihood is 
|	recomputed. If the edge of `focal_node' has been changed, the function invalidateBothEnds(`focal_node') should also
|	be called to invalidate the one remaining CLA not invalidated by this function.
*/
void TreeLikelihood::invalidateAwayFromNode(
  TreeNode & focal_node)		/**< invalidation of conditional likelihood arrays will proceed outward from this node */
	{
	NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::invalidateNode, this, _1, _2);
	effective_postorder_edge_iterator(&focal_node, validFunctor); // need only construct unnamed iterator object
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Unconditionally discards cached parental (and, if ref_nd is internal, filial) conditional likelihood arrays if
|	present. This function always returns false so it can be used in conjunction with effective_postorder_edge_iterator
|	to discard every cached CLA in the entire tree (something that needs to be done when any move is accepted).
*/
bool TreeLikelihood::discardCacheBothEnds(TreeNode * ref_nd, TreeNode * /* unused */)
	{
	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();

		// Remove cached parental CLAs if they exist
		if (td->parCachedCLA)
			{
			cla_pool.putCondLikelihood(td->parCachedCLA);
			td->parCachedCLA.reset();
			}
		}
	else
		{
		// Internal nodes have both parental and filial CLAs
		InternalData * id = ref_nd->GetInternalData();

		// Remove cached parental CLAs if they exist
		if (id->parCachedCLA)
			{
			cla_pool.putCondLikelihood(id->parCachedCLA);
			id->parCachedCLA.reset();
			}

		// Remove cached filial CLAs if they exist
		if (id->childCachedCLA)
			{
			cla_pool.putCondLikelihood(id->childCachedCLA);
			id->childCachedCLA.reset();
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores parental conditional likelihood arrays from cache and stores the existing conditional likelihood arrays 
|	for later use. This function always returns false so it can be used in conjunction with 
|	effective_postorder_edge_iterator to restore every parental CLA from cache in the entire tree. If `ref_nd' 
|	has a working parental CLA but no corresponding cached parental CLA, then the working parental CLA is discarded. 
|	This is done because this working parental CLA must have been calculated while the tree was in a proposed state, 
|	and thus would be invalid if the tree were reverted.
*/
bool TreeLikelihood::restoreFromCacheParentalOnly(TreeNode * ref_nd, TreeNode * /* unused */)
	{
	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();

		// Store working CLA in any case
		if (td->parWorkingCLA)
			cla_pool.putCondLikelihood(td->parWorkingCLA);
		td->parWorkingCLA.reset();

		// Move cached to working if cached exists
		if (td->parCachedCLA)
			{
			td->parWorkingCLA = td->parCachedCLA;
			td->parCachedCLA.reset();
			}
		}
	else
		{
		InternalData * id = ref_nd->GetInternalData();

		// Store working CLA in any case
		if (id->parWorkingCLA)
			cla_pool.putCondLikelihood(id->parWorkingCLA);
		id->parWorkingCLA.reset();

		// Move cached to working if cached exists
		if (id->parCachedCLA)
			{
			id->parWorkingCLA = id->parCachedCLA;
			id->parCachedCLA.reset();
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores parental (and, if ref_nd is internal, filial) conditional likelihood arrays from cache and stores the
|	existing conditional likelihood arrays for later use. This function always returns false so it can be used in 
|	conjunction with effective_postorder_edge_iterator to restore every CLA from cache in the entire tree. If `ref_nd' 
|	has a working CLA but no corresponding cached CLA, then the working CLA is discarded. This is done because this 
|	working CLA must have been calculated while the tree was in a proposed state, and thus would be invalid if the tree
|	were reverted.
*/
bool TreeLikelihood::restoreFromCacheBothEnds(TreeNode * ref_nd, TreeNode * /* unused */)
	{
	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();

		// Store working CLA in any case
		if (td->parWorkingCLA)
			cla_pool.putCondLikelihood(td->parWorkingCLA);
		td->parWorkingCLA.reset();

		// Move cached to working if cached exists
		if (td->parCachedCLA)
			{
			td->parWorkingCLA = td->parCachedCLA;
			td->parCachedCLA.reset();
			}
		}
	else
		{
		// Internal nodes have both parental and filial CLAs
		InternalData * id = ref_nd->GetInternalData();

		// Restore parental CLA from cache

		// Store working CLA in any case
		if (id->parWorkingCLA)
			cla_pool.putCondLikelihood(id->parWorkingCLA);
		id->parWorkingCLA.reset();

		// Move cached to working if cached exists
		if (id->parCachedCLA)
			{
			id->parWorkingCLA = id->parCachedCLA;
			id->parCachedCLA.reset();
			}

		// Restore filial CLA from cache

		// Store working CLA in any case
		if (id->childWorkingCLA)
			cla_pool.putCondLikelihood(id->childWorkingCLA);
		id->childWorkingCLA.reset();

		// Move cached to working if cached exists
		if (id->childCachedCLA)
			{
			id->childWorkingCLA = id->childCachedCLA;
			id->childCachedCLA.reset();
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function used in conjunction with effective_postorder_edge_iterator to restore the conditional likelihood arrays 
|	(CLAs) of all nodes starting with the supplied focal node `ref_nd'. For nodes that are descendants of the focal node
|	(descendant means that a node can be found using only left child and right sib pointers starting from the focal
|	node), it is the parental CLAs that are restored. For a node that is an ancestor of the focal node, it is the 
|	filial CLAs that are restored. For all other nodes (e.g. on independent lineages derived from an ancestral node), 
|	it is the parental CLAs that are restored. This function always returns false because the goal is to restore every 
|	node that needs to be invalidated. If a node is has a working CLA but no corresponding cached CLA, then the working
|	CLA is discarded. This is done because this working CLA must have been calculated while the tree was in a proposed
|	state, and thus would be invalid if the tree were reverted.
*/
bool TreeLikelihood::restoreFromCacheNode(TreeNode * ref_nd, TreeNode * neighbor_closer_to_likelihood_root)
	{
	if (ref_nd->IsTip() && !ref_nd->IsTipRoot())
		{
		TipData * td = ref_nd->GetTipData();

		// Store the working CLA in any case
		if (td->parWorkingCLA)
			cla_pool.putCondLikelihood(td->parWorkingCLA);
		td->parWorkingCLA.reset();

		// Move cached to working if there is a cached CLA
		if (td->parCachedCLA)
			{
			td->parWorkingCLA = td->parCachedCLA;
			td->parCachedCLA.reset();
			}
		}
	else
		{
		if (ref_nd->GetParent() == neighbor_closer_to_likelihood_root)
			{
			InternalData * id = ref_nd->GetInternalData();

			// Store the working CLA in any case
			if (id->parWorkingCLA)
				cla_pool.putCondLikelihood(id->parWorkingCLA);
			id->parWorkingCLA.reset();

			// Move cached to working if there is a cached CLA
			if (id->parCachedCLA)
				{
				id->parWorkingCLA = id->parCachedCLA;
				id->parCachedCLA.reset();
				}
			}
		else
			{
			// neighbor_closer_to_likelihood_root is the actual child
			if (neighbor_closer_to_likelihood_root->IsTip())
				return false;

			// Either ref_nd is not a tip, or it is the tip serving as the root node
			InternalData * id = neighbor_closer_to_likelihood_root->GetInternalData();

			// Store the working CLA in any case
			if (id->childWorkingCLA)
				cla_pool.putCondLikelihood(id->childWorkingCLA);
			id->childWorkingCLA.reset();

			// Move cached to working if there is a cached CLA
			if (id->childCachedCLA)
				{
				id->childWorkingCLA = id->childCachedCLA;
				id->childCachedCLA.reset();
				}
			}
		}
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores all conditional likelihood arrays (CLAs) pointing toward the supplied `focal_node' from cache. That is, the
|	parental CLA of a child of `focal_node' would be restored, as would be the filial CLA for the parent of the 
|	`focal_node'. It may be necessary to call the function invalidateBothEnds(`focal_node') also to restore the one 
|	remaining CLA not restored by this function.
*/
void TreeLikelihood::restoreFromCacheAwayFromNode(
  TreeNode & focal_node)		/**< restoration of cached conditional likelihood arrays will proceed outward from this node */
	{
	NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::restoreFromCacheNode, this, _1, _2);
	effective_postorder_edge_iterator(&focal_node, validFunctor); // need only construct unnamed iterator object
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Discards all cached conditional likelihood arrays (CLAs) pointing toward the supplied `focal_node'. That is, the
|	cached parental CLA of a child of `focal_node' would be discarded, as would be the cached filial CLA for the parent 
|	of the `focal_node'. None of the working CLAs are affected. It may be necessary to call the function 
|	discardCacheBothEnds(`focal_node') also to discard the cache from the two remaining CLAs not affected by this 
|	function.
*/
void TreeLikelihood::discardCacheAwayFromNode(
  TreeNode & focal_node)		/**< restoration of cached conditional likelihood arrays will proceed outward from this node */
	{
	NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::discardCacheBothEnds, this, _1, _2);
	effective_postorder_edge_iterator(&focal_node, validFunctor); // need only construct unnamed iterator object
	}

#if 0
class CalcTransitionMatrixForOneNode : public std::unary_function<TreeNode &, void>
	{
	private:
		TreeLikelihood & treelike;

	public:
		CalcTransitionMatrixForOneNode(TreeLikelihood & tl) : treelike(tl) {}
		void operator()(TreeNode & nd)
			{
#			error do not use unless root node special case is taken into account
			if (!nd.IsTipRoot())
				{
				double edge_len = nd.GetEdgeLen();
				if (nd.IsTip())
					{
					TipData & ndTD = *(nd.GetTipData());
					treelike.calcTMatForSim(ndTD, edge_len);
					}
				else
					{
					InternalData & ndID	= *(nd.GetInternalData());
					treelike.calcPMat(ndID, edge_len);
					}
				}
			}
	};
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Simulates data for `nchar' characters using the current model and edge lengths. Transition matrices are recomputed
|	if `refresh_probs' is true; if `refresh_probs' is false, assumes transition probability matrices are up-to-date. 
|	Only the `num_states' primary states are generated, and the state created for each node is stored initially in the 
|	`state' data member of the TipData or InternalData structure associated with the tip or internal node, respectively.
|	This function expects these data structures to be in place (the TreeLikelihood::prepareForSimulation and 
|	TreeLikelihood::prepareForLikelihood functions both perform this task). As states are generated for tip nodes, they
|	are copied into a temporary pattern inside `sim_data', and this pattern is then inserted into a pattern map inside 
|	`sim_data' when it is completed. Assumes that `nchar' is greater than zero and that the shared pointers `sim_data' 
|	and `rng' actually point to real objects.
|	
|	The following makes use of T matrices, which are transposed, augmented transition matrices. These are transposed 
|	because the "from" states form the columns rather than the rows. They are augmented because, ordinarily, there are
|	additional rows corresponding to ambiguities observed in some tip nodes. With simulated data, however, there are 
|	never any ambiguities, so in this case T matrices are nothing more than transposed transition matrices. 
|	Important: note that T matrices are only used in the TipData structures; InternalData structures store normal 
|	untransposed transition probability matrices in which the rows form the "from" states. 
|	
|	Here is an example of a T matrix created using the JC model and an edge length equal to 0.1:
|>	
|	            |--------------- from state -------------|
|	                0          1          2          3
|	  	            A          C          G          T
|	t  0   A     0.90638    0.03121    0.03121    0.03121
|	o  1   C     0.03121    0.90638    0.03121    0.03121
|	   2   G     0.03121    0.03121    0.90638    0.03121
|	s  3   T     0.03121    0.03121    0.03121    0.90638
|	t  4   N     1.00000    1.00000    1.00000    1.00000 \
|	a  5 {GT}    0.06241    0.06241    0.93759    0.93759  | These rows not present if prepareForSimulation function
|	t  6 {ACT}   0.96879    0.96879    0.09362    0.96879  | was used to create the TipData structures
|	e  7 {AG}    0.93757    0.06241    0.93759    0.06241 /
|>
|	The `pMatrixTranspose' data member in TipData structures holds the array of T matrices (one T matrix for each 
|	rate category).
*/
void TreeLikelihood::simulateImpl(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar, bool refresh_probs)
	{
	PHYCAS_ASSERT(sim_data);
	PHYCAS_ASSERT(rng);
	PHYCAS_ASSERT(nchar > 0);

	// Recalculate transition probabilities if requested
	//@POL using for_each would simplify this
	if (refresh_probs)
		{
		preorder_iterator nd = t->begin();

		// First preorder node is the root node and represents a special case
		// Its transition matrices must be computed using the "subroot" node's edge length
		// The subroot node's transition matrices need not be calculated 
		TipData & ndTD = *(nd->GetTipData());
		TreeNode * subroot = nd->GetLeftChild();
		PHYCAS_ASSERT(subroot);
		PHYCAS_ASSERT(!subroot->GetRightSib());	//@POL need to create a IsSubroot() member function for TreeNode
		calcTMatForSim(ndTD, subroot->GetEdgeLen());
		++nd;

		// Skip subroot node as its transition matrices are never used and thus do not need to be computed
		++nd;

		// Process the remaining nodes in the tree
		for (; nd != t->end(); ++nd)
			{
			if (nd->IsTip())
				{
				TipData & ndTD = *(nd->GetTipData());
				calcTMatForSim(ndTD, nd->GetEdgeLen());
				}
			else
				{
				InternalData & ndID	= *(nd->GetInternalData());
				calcPMat(ndID.getPMatrices(), nd->GetEdgeLen());
				}
			}
		}

	// Create a vector of cumulative state frequencies to use in choosing starting states
	const std::vector<double> & freqs = model->getStateFreqs();
	std::vector<double> cum_freqs(num_states, 0.0);
	std::partial_sum(freqs.begin(), freqs.end(), cum_freqs.begin());

	// Create a vector of cumulative rate probabilities to use in choosing relative rates
	std::vector<double> cum_rate_probs(num_rates, 0.0);
	std::partial_sum(rate_probs.begin(), rate_probs.end(), cum_rate_probs.begin());

	sim_data->resetPatternLength(t->GetNTips());
	sim_data->wipePattern();

#if 0  //POL temporary debugging section BEGIN
	//std::ofstream doof("doof_check.txt", std::ios::out | std::ios::app); 
	std::ofstream doof("doof_check.txt"); 
	doof << "\nNEW DATA SET BEGINNING" << std::endl;
	doof << "model parameter names:  " << model->paramHeader() << std::endl;
	doof << "model parameter values: " << model->paramReport() << std::endl;

	// Go through all nodes and show their edge lengths and transition probability matrices
	{
		preorder_iterator nd = t->begin();

		// First preorder node is the root node and represents a special case
		// Its transition matrices must be computed using the "subroot" node's edge length
		// The subroot node's transition matrices need not be calculated 
		TipData & ndTD = *(nd->GetTipData());
		double * * * p = ndTD.getTransposedPMatrices();
		TreeNode * subroot = nd->GetLeftChild();
		doof << "Subroot node edge length = " << subroot->GetEdgeLen() << std::endl;
		doof << "Transposed transition probability matrices:" << std::endl;
		for (unsigned rr = 0; rr < num_rates; ++rr)
			{
			doof << "  Relative rate " << rr << " = " << rate_means[rr] << std::endl;
			doof << str(boost::format("  %12.5f %12.5f %12.5f %12.5f\n") % p[rr][0][0] % p[rr][0][1] % p[rr][0][2] % p[rr][0][3]);
			doof << str(boost::format("  %12.5f %12.5f %12.5f %12.5f\n") % p[rr][1][0] % p[rr][1][1] % p[rr][1][2] % p[rr][1][3]);
			doof << str(boost::format("  %12.5f %12.5f %12.5f %12.5f\n") % p[rr][2][0] % p[rr][2][1] % p[rr][2][2] % p[rr][2][3]);
			doof << str(boost::format("  %12.5f %12.5f %12.5f %12.5f\n") % p[rr][3][0] % p[rr][3][1] % p[rr][3][2] % p[rr][3][3]);
			doof << std::endl;
			}
		++nd;

		// Skip subroot node as its transition matrices are never used and thus do not need to be computed
		++nd;

		// Process the remaining nodes in the tree
		for (; nd != t->end(); ++nd)
			{
			if (nd->IsTip())
				{
				TipData & ndTD = *(nd->GetTipData());
				double * * * p = ndTD.getTransposedPMatrices();
				doof << "Tip node " << nd->GetNodeNumber() << " edge length = " << nd->GetEdgeLen() << std::endl;
				doof << "Transposed transition probability matrices:" << std::endl;
				for (unsigned rr = 0; rr < num_rates; ++rr)
					{
					doof << "  Relative rate " << rr << " = " << rate_means[rr] << std::endl;
					doof << str(boost::format("  %12.5f %12.5f %12.5f %12.5f\n") % p[rr][0][0] % p[rr][0][1] % p[rr][0][2] % p[rr][0][3]);
					doof << str(boost::format("  %12.5f %12.5f %12.5f %12.5f\n") % p[rr][1][0] % p[rr][1][1] % p[rr][1][2] % p[rr][1][3]);
					doof << str(boost::format("  %12.5f %12.5f %12.5f %12.5f\n") % p[rr][2][0] % p[rr][2][1] % p[rr][2][2] % p[rr][2][3]);
					doof << str(boost::format("  %12.5f %12.5f %12.5f %12.5f\n") % p[rr][3][0] % p[rr][3][1] % p[rr][3][2] % p[rr][3][3]);
					doof << std::endl;
					}
				}
			else
				{
				InternalData & ndID	= *(nd->GetInternalData());
				double * * * p = ndID.getPMatrices();
				doof << "Internal node " << nd->GetNodeNumber() << " edge length = " << nd->GetEdgeLen() << std::endl;
				doof << "Transition probability matrices:" << std::endl;
				for (unsigned rr = 0; rr < num_rates; ++rr)
					{
					doof << "  Relative rate " << rr << " = " << rate_means[rr] << std::endl;
					doof << str(boost::format("  %12.5f %12.5f %12.5f %12.5f\n") % p[rr][0][0] % p[rr][0][1] % p[rr][0][2] % p[rr][0][3]);
					doof << str(boost::format("  %12.5f %12.5f %12.5f %12.5f\n") % p[rr][1][0] % p[rr][1][1] % p[rr][1][2] % p[rr][1][3]);
					doof << str(boost::format("  %12.5f %12.5f %12.5f %12.5f\n") % p[rr][2][0] % p[rr][2][1] % p[rr][2][2] % p[rr][2][3]);
					doof << str(boost::format("  %12.5f %12.5f %12.5f %12.5f\n") % p[rr][3][0] % p[rr][3][1] % p[rr][3][2] % p[rr][3][3]);
					doof << std::endl;
					}
				}
			}
	}
	doof.close();
#endif	//POL temporary debugging section END

	for (unsigned character = 0; character < nchar; ++character)
		{
		// Choose a rate for this character (actually, choose index, the actual rate is rate_means[r])
		unsigned r = 0;
		if (num_rates > 1)
			{
			// warning: removing the if statement will invalidate all examples involving simulated data with rate
			// homogeneity because of the call to rng->Uniform here!
			r = (unsigned)(std::lower_bound(cum_rate_probs.begin(), cum_rate_probs.end(), rng->Uniform()) - cum_rate_probs.begin());
			}

		// Generate the starting state
		int8_t j = (unsigned)(std::lower_bound(cum_freqs.begin(), cum_freqs.end(), rng->Uniform()) - cum_freqs.begin());

		// Assign starting state to the tip node currently serving as the root of the tree
		preorder_iterator nd = t->begin();
		TipData & rootTD = *(nd->GetTipData());
		rootTD.state = j;

		sim_data->setState(nd->GetNodeNumber(), j);

		// Go ahead and generate the state for the (only) descendant of the root node (the "subroot" node)
		// Note that the root node's T matrix is used for this calculation; the P matrix of the subroot node
		// is never computed
		unsigned parent_state = (unsigned)j;

		// Get the T matrix for the tip node serving as the root
		double * * Tmatrix = rootTD.pMatrixTranspose[r];

		// Choose a uniform random deviate
		double u = rng->Uniform();

		// Spin the roulette wheel to choose a state for the subroot node
		double cum = 0.0;
		unsigned i = 0;
		for (; i < num_states; ++i)
			{
			double pr = Tmatrix[i][parent_state];
			//std::cerr << str(boost::format("Tmatrix[%d][%d] = %f") % i % parent_state % pr) << std::endl;
			cum += pr;
			if (u < cum)
				break;
			}

		// Increment iterator so that nd now refers to the subroot (sole descendant of the root)
		++nd;

		// Assign the new state to the subroot node
		InternalData & ndID = *(nd->GetInternalData());
		ndID.state = (int8_t)i;
		//std::cerr << "  Assigning state " << i << " to node " << nd->GetNodeNumber() << std::endl;

		// Walk the remainder of the tree using the preorder sequence, generating data for each node along the way
		for (++nd; nd != t->end(); ++nd)
			{
			// Get state of parent of nd
			TreeNode * parent = nd->GetParent();
			parent_state = UINT_MAX;
			if (parent->IsTip())
				{
				TipData * parentTD = parent->GetTipData();
				parent_state = (unsigned)parentTD->state;
				}
			else
				{
				InternalData * parentID = parent->GetInternalData();
				parent_state = (unsigned)parentID->state;
				}
			PHYCAS_ASSERT(parent_state < num_states);

			if (nd->IsTip())
				{
				// Get the T matrix
				TipData & ndTD = *(nd->GetTipData());
				double * * Tmatrix = ndTD.pMatrixTranspose[r];

				// Choose a uniform random deviate
				double u = rng->Uniform();

				// Spin the roulette wheel and assign a state to nd
				double cum = 0.0;
				unsigned i = 0;
				for (; i < num_states; ++i)
					{
					double pr = Tmatrix[i][parent_state];
					//std::cerr << str(boost::format("Tmatrix[%d][%d] = %f") % i % parent_state % pr) << std::endl;
					cum += pr;
					if (u < cum)
						break;
					}
				ndTD.state = (int8_t)i;
				sim_data->setState(nd->GetNodeNumber(), (int8_t)i);

				//std::cerr << "  Assigning state " << i;
				}
			else
				{
				// Get the T matrix
				InternalData & ndID = *(nd->GetInternalData());
				double * * Pmatrix = ndID.pMatrices[r];

				// Choose a uniform random deviate
				double u = rng->Uniform();

				// Spin the roulette wheel and assign a state to nd
				double cum = 0.0;
				unsigned i = 0;
				for (; i < num_states; ++i)
					{
					double pr = Pmatrix[parent_state][i];
					//std::cerr << str(boost::format("Pmatrix[%d][%d] = %f") % parent_state % i % pr) << std::endl;
					cum += pr;
					if (u < cum)
						break;
					}
				ndID.state = (int8_t)i;
				//std::cerr << "  Assigning state " << i;
				}

			//std::cerr << " to node " << nd->GetNodeNumber() << std::endl;
			}
		//std::cerr << std::endl;

		// We are now finished simulating data for one character, so insert the pattern just generated
		// into the pattern map maintained by sim_data; the 1.0 means that the count for this pattern
		// should be incremented by 1
		sim_data->insertPattern(1.0);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Stores all CondLikelihood objects stored in TipData or InternalData structures. This will require all CLAs to be
|	recomputed the next time the likelihood needs to be computed. Returns the subroot node (only child of root node).
*/
TreeNode * TreeLikelihood::storeAllCLAs(
  TreeShPtr t)				/**< is the tree from which all CLAs are to be removed */
	{
	// Start at the root node
	TreeNode * nd = t->GetFirstPreorder();
	PHYCAS_ASSERT(nd);

	// Move to the subroot node
	nd = nd->GetNextPreorder();
	PHYCAS_ASSERT(nd);

	// Invalidate (and do not cache) all CLAs from the tree. This will require all CLAs to be recomputed 
	// when the likelihood is computed using the subroot as the likelihood root. This path should be taken
	// if a parameter is changed that invalidates the entire tree.
	NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::invalidateBothEndsDiscardCache, this, _1, _2);
	effective_postorder_edge_iterator(nd, validFunctor); // constructor does all the work we need
	invalidateBothEndsDiscardCache(nd);

	return nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sanity check to make sure TreeLikelihood::storeAllCLAs really has stored all CLAs. Returns true if any CLAs are 
|	found, false if no CLAs are found (not even in cache positions).
*/
bool TreeLikelihood::debugCheckCLAsRemainInTree(
  TreeShPtr t) const	/**< is the tree to check */
	{
	for (TreeNode * nd = t->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			{
			TipData * td = nd->GetTipData();
			if (td != NULL)
				{
				if (td->parWorkingCLA)
					return true;
				if (td->parCachedCLA)
					return true;
				}
			}
		else
			{
			InternalData * id = nd->GetInternalData();
			if (id != NULL)
				{
				if (id->parWorkingCLA)
					return true;
				if (id->parCachedCLA)
					return true;

				if (id->childWorkingCLA)
					return true;
				if (id->childCachedCLA)
					return true;
				}
			}
		}	// preorder loop over all nodes

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This is the function that needs to be called to recompute the log-likelihood. If `likelihood_root' is not NULL, 
|	the calculation will use that node as the root of the likelihood calculation, and it will be assumed that all
|	conditional likelihood arrays (CLAs) are correctly calculated or have been invalidated if they need to be 
|	recomputed. On the other hand, if `likelihood_root' is NULL, then all CLAs will be invalidated and recomputed
|	(useful if a parameter has been changed that requires recalculation of all CLAs).
*/
double TreeLikelihood::calcLnL(
  TreeShPtr t)
	{
	if (no_data)
		return 0.0;

	// Compute likelihood using likelihood_root if specified
	// Assume that if likelihood_root has been specified, then the necessary 
	// CLA invalidations have already been performed.
	TreeNode * nd = likelihood_root;
	if (nd == NULL)
		{
#if 0
		// If no likelihood_root has been specified, use the subroot node (and
		// invalidate the entire tree to be safe)
		nd = t->GetFirstPreorder();
		PHYCAS_ASSERT(nd);

		// Move to the subroot node
		nd = nd->GetNextPreorder();
		PHYCAS_ASSERT(nd);

		// The subroot node is the new likelihood_root
		likelihood_root = nd;

		// Invalidate (and do not cache) all CLAs from the tree. This will require all CLAs to be recomputed 
		// when the likelihood is computed using the subroot as the likelihood root. This path should be taken
		// if a parameter is changed that invalidates the entire tree.
		NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::invalidateBothEndsDiscardCache, this, _1, _2);
		effective_postorder_edge_iterator(nd, validFunctor); // constructor does all the work we need
		invalidateBothEndsDiscardCache(nd);
#else
		// If no likelihood_root has been specified, invalidate the entire tree to be safe
		nd = storeAllCLAs(t);

		// The subroot node will be the new likelihood_root
		likelihood_root = nd;
#endif
		}

	PHYCAS_ASSERT(nd);
	PHYCAS_ASSERT(nd->IsInternal());

	// Calculate log-likelihood using nd as the likelihood root
	double lnL = calcLnLFromNode(*nd);

	return lnL;
	}

void TreeLikelihood::debugSaveCLAs(TreeShPtr t, std::string fn, bool overwrite)
	{
	std::ofstream tmpf;
	if (overwrite)
		tmpf.open(fn.c_str());
	else
		{
		tmpf.open(fn.c_str(), std::ios::out | std::ios::app);
		}

	for (TreeNode * nd = t->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			{
			tmpf << "\n\nTip node number " << nd->GetNodeNumber();
			if (nd->IsTipRoot())
                tmpf << "\n  parent = None";
			else
                tmpf << "\n  parent = " << nd->GetParent()->GetNodeNumber();
            tmpf << "\n  children = ";
			for (TreeNode * child = nd->GetLeftChild(); child != NULL; child = child->GetRightSib())
                tmpf << child->GetNodeNumber() << " ";

			TipData * td = nd->GetTipData();
			if (td->parentalCLAValid())
				{
				CondLikelihoodShPtr cla = td->getParentalCondLikePtr();
				unsigned sz = cla->getCLASize();
				LikeFltType * arr = cla->getCLA();
	            tmpf << "\n  cla length = " << sz;
	            tmpf << "\n  parental = ";
				for (unsigned i = 0; i < sz; ++i)
					{
					tmpf << str(boost::format("%20.12f ") % arr[i]);
					}
				}
			else
				{
	            tmpf << "\n  parental CLA = None";
				}
			}
		else
			{
			tmpf << "\n\nInternal node number " << nd->GetNodeNumber();

            tmpf << "\n  parent = " << nd->GetParent()->GetNodeNumber();
            tmpf << "\n  children = ";
			for (TreeNode * child = nd->GetLeftChild(); child != NULL; child = child->GetRightSib())
                tmpf << child->GetNodeNumber() << " ";

			InternalData * id = nd->GetInternalData();
			if (id->parentalCLAValid())
				{
				CondLikelihoodShPtr cla = id->getParentalCondLikePtr();
				unsigned sz = cla->getCLASize();
				LikeFltType * arr = cla->getCLA();
	            tmpf << "\n  parental CLA length = " << sz;
	            tmpf << "\n  parental CLA = ";
				for (unsigned i = 0; i < sz; ++i)
					{
					tmpf << str(boost::format("%20.12f ") % arr[i]);
					}
				}
			else
				{
	            tmpf << "\n  parental CLA = None";
				}

			if (id->filialCLAValid())
				{
				CondLikelihoodShPtr cla = id->getChildCondLikePtr();
				unsigned sz = cla->getCLASize();
				LikeFltType * arr = cla->getCLA();
	            tmpf << "\n  filial CLA length = " << sz;
	            tmpf << "\n  filial CLA = ";
				for (unsigned i = 0; i < sz; ++i)
					{
					tmpf << str(boost::format("%20.12f ") % arr[i]);
					}
				}
			else
				{
	            tmpf << "\n  filial CLA = None";
				}
			}
		}

	tmpf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the log-likelihood using `focal_node' as the likelihood root (i.e. the root of the likelihood calculation,
|	but not necessarily the root of the tree).
*/
double TreeLikelihood::calcLnLFromNode(
  TreeNode & focal_node)	/**< is the likelihood root (i.e. node around which the likelihood will be computed) */
	{
	double lnL;
	if (no_data)
		lnL =  0.0;
	else
		{
		PHYCAS_ASSERT(!focal_node.IsTip());

		// valid_functor will return true if the conditional likelihood arrays pointing away from the
		// focal_node are up-to-date, false if they need to be recomputed
		NodeValidityChecker valid_functor = boost::bind(&TreeLikelihood::isValid, this, _1, _2);

		// iter will visit nodes that need their CLAs updated centripetally (like a postorder traversal 
		// but also coming from below the focal node). Each node visited is guaranteed by valid_functor 
		// to need its CLA updated.

		effective_postorder_edge_iterator iter(&focal_node, valid_functor);
		effective_postorder_edge_iterator iter_end;
		for (; iter != iter_end; ++iter)
			{
			refreshCLA(*iter->first, iter->second);
			}

		// We have now brought all neighboring CLAs up-to-date, so we can now call harvestLnL to
		// compute the likelihood
		EdgeEndpoints edge(&focal_node, NULL);
		lnL = harvestLnL(edge);
		}
	++nevals;
	return lnL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates the TipData data structure needed to store the data for one tip (the tip corresponding to the supplied
|	`row' index in the data matrix `mat'). Returns a pointer to the newly-created TipData structure. See documentation 
|	for the TipData structure for more explanation.
*/
TipData * TreeLikelihood::allocateTipData(  //POLBM TreeLikelihood::allocateTipData
  unsigned row) 	/**< is the row of the data matrix corresponding to the data for this tip node */
	{
	std::map<int8_t, int8_t>					globalToLocal;
	std::vector<unsigned int>					stateListVec;
	std::map<int8_t, int8_t>::const_iterator	foundElement;

	int8_t *									tipSpecificStateCode	= new int8_t[(unsigned)pattern_map.size()];
	//@POL 21-Nov-2005 make tipSpecificStateCode a shared_array or a std:Vector - currently I don't think these are being deleted

	const int8_t								ns						= num_states;
	const int8_t								nsPlusOne				= num_states + 1;
	unsigned									nPartialAmbig			= 0;

	// Loop through all patterns for this row of the matrix. For each global state code encountered,
	// determine which local state code it represents, and build up the tipSpecificStateCode array
	// as we go.
	unsigned i = 0;
#if 0
	if (model->isCodonModel())
		{
		// Read three nucleotides at a time and interpret the triplet as a codon state
		//@POL Currently, any ambiguity at any codon position results in a missing data entry for the codon state
		for (PatternMapType::const_iterator it = pattern_map.begin(); it != pattern_map.end(); ++it, ++i)
			{
			const int8_t globalStateCode = (it->first)[row];

			if (globalStateCode < nsPlusOne)
				{
				// no partial ambiguity, but may be gap state
				tipSpecificStateCode[i] = (globalStateCode < 0 ? ns : globalStateCode);
				}
			else
				{
				// partial ambiguity
				foundElement = globalToLocal.find(globalStateCode);
				if (foundElement == globalToLocal.end())
					{
					// state code needs to be added to map
					globalToLocal[globalStateCode] = nPartialAmbig + nsPlusOne;
					stateListVec.push_back(state_list_pos[globalStateCode]);
					tipSpecificStateCode[i] = nPartialAmbig + nsPlusOne;
					nPartialAmbig++;
					}
				else
					{
					// state code is already in the map
					tipSpecificStateCode[i] = foundElement->second;
					}
				}
			}
		}
	else
		{
#endif
		for (PatternMapType::const_iterator it = pattern_map.begin(); it != pattern_map.end(); ++it, ++i)
			{
			const int8_t globalStateCode = (it->first)[row];

			if (globalStateCode < nsPlusOne)
				{
				// no partial ambiguity, but may be gap state
				tipSpecificStateCode[i] = (globalStateCode < 0 ? ns : globalStateCode);
				}
			else
				{
				// partial ambiguity
				foundElement = globalToLocal.find(globalStateCode);
				if (foundElement == globalToLocal.end())
					{
					// state code needs to be added to map
					globalToLocal[globalStateCode] = nPartialAmbig + nsPlusOne;
					stateListVec.push_back(state_list_pos[globalStateCode]);
					tipSpecificStateCode[i] = nPartialAmbig + nsPlusOne;
					nPartialAmbig++;
					}
				else
					{
					// state code is already in the map
					tipSpecificStateCode[i] = foundElement->second;
					}
				}
			}
#if 0
		}
#endif

#	if 0
		//POL-debug
		std::map<int8_t, int8_t>::const_iterator mapiter; 
		unsigned z;

		//std::cerr << "\n\nInside allocateTipData function:" << std::endl;

		//std::cerr << "\n  ns: " << (int)ns << std::endl;

		std::cerr << "\n  Adding TipData for node number " << row << ": ";
				for (z = 0; z < num_patterns; ++z)
					std::cerr << ' ' << (int)(myRow[z]);
		std::cerr << std::endl;

		//std::cerr << "\n  Contents of state_list_pos:" << std::endl;
		//for (vector<unsigned int>::const_iterator ziter = state_list_pos.begin(); ziter != state_list_pos.end(); ++ziter)
		//	std::cerr << "\n    " << (*ziter);
		//std::cerr << std::endl;

		//std::cerr << "\n  nPartialAmbig: " << nPartialAmbig << std::endl;

		//std::cerr << "\n  Contents of tipSpecificStateCode:" << std::endl;
		//for (z = 0; z < nPatterns; ++z)
		//	std::cerr << "\n    " << (int)(tipSpecificStateCode[z]);
		//std::cerr << std::endl;

		//if (globalToLocal.empty())
		//	std::cerr << "\n  globalToLocal map is empty" << std::endl;
		//else
		//	{
		//	std::cerr << "\n  Contents of globalToLocal map:" << std::endl;
		//	for (mapiter = globalToLocal.begin(); mapiter != globalToLocal.end(); ++mapiter)
		//		std::cerr << "\n    " << (int)(mapiter->first) << " (global) -> " << (int)(mapiter->second) << " (local)";
		//	std::cerr << std::endl;
		//	}

#	endif

	return new TipData(	stateListVec,												// stateListPosVec
						boost::shared_array<const int8_t>(tipSpecificStateCode),	// stateCodesShPtr
						num_rates,													// number of relative rate categories
						num_states,													// number of states in the model
						NULL,														// pMatTranspose
						true,														// managePMatrices
						cla_pool);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates the InternalData data structure needed to store the conditional likelihood arrays and transition matrices
|	needed for likelihood calcuations. Returns a pointer to a newly-constructed InternalData structure. See the 
|	documentation for InternalData for more explanation.
*/
InternalData * TreeLikelihood::allocateInternalData()
	{
	return new InternalData(num_patterns,				// number of site patterns
							num_rates,					// number of relative rate categories
							num_states,					// number of model states
							NULL,						// pMat
							true,						// managePMatrices
							cla_pool);					
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function that deletes the TipData structure allocated by the allocateTipData() function. The intention is for this
|	function to be given to TreeNode objects (in the form of a boost::function object) whenever allocateTipData() is 
|	used to allocate memory for a TipData structure. The TreeNode destructor can then use the boost::function object to
|	delete the memory required by the TipData structure without knowing any details of the TipData structure (i.e.,
|	the file that defines the TreeNode destructor does not need to include a header file that declares the details of
|	the TipData class.
*/
void deallocateTipData(TipData * p)
	{
	delete p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function that deletes the InternalData structure allocated by the allocateInternalData() function. The intention is 
|	for this function to be given to TreeNode objects (in the form of a boost::function object) whenever 
|	allocateInternalData() is used to allocate memory for an InternalData structure. The TreeNode destructor can then 
|	use the boost::function object to delete the memory required by the InternalData structure without knowing any 
|	details of the InternalData structure (i.e., the file that defines the TreeNode destructor does not need to include 
|	a header file that declares the details of the InternalData class.
*/
void deallocateInternalData(InternalData * p)
	{
	delete p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a basic TipData object for all tip nodes and calls allocateInternalData() for all internal nodes in the 
|	supplied Tree. This function does not require a data matrix and does not add observed data to TipData structures.
|	Silently returns if root node already has a TipData structure. This is taken to indicate that prepareForLikelihood
|	or prepareForSimulation was previously called for this tree, in which case all data structures needed for 
|	simulation are already present.
*/
void TreeLikelihood::prepareForSimulation(
  TreeShPtr t)			/**< is the tree to decorate */
	{
	TreeNode::TipDataDeleter		td_deleter	= &deallocateTipData;
	TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;

	preorder_iterator nd = t->begin();

	// Only proceed if root node does not already have a TipData structure
	if (nd->GetTipData() != NULL)
		return;

	for (; nd != t->end(); ++nd)
		{
		if (nd->IsTip())
			{
			TipData * td = 	new TipData(num_rates,	num_states, cla_pool);	//@POL should be using shared_ptr here?
			nd->SetTipData(td, td_deleter);
			}
		else
			{
			InternalData * cl = allocateInternalData();
			nd->SetInternalData(cl, cl_deleter);
			}
		}
	}

#if defined(INTERFACE_WITH_CIPRES)
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void TreeLikelihood::copyDataFromIDLMatrix(
  const PhycasIDLishMatrix & m)	/**< */
	{
	// Copy information about state codes 
	state_list.clear();
	state_list_pos.clear();
	int curr_pos = 0;
	typedef std::vector<CIPR_State_t> CIPRStateVect;
	std::vector<CIPRStateVect>::const_iterator rowit = m.state_lookup.begin();
	for (unsigned i = 0; rowit != m.state_lookup.end(); ++rowit, ++i)
		{
		state_list_pos.push_back(curr_pos);
		curr_pos += (unsigned)rowit->size();
		std::vector<CIPR_State_t>::const_iterator it = rowit->begin();
		for (; it != rowit->end(); ++it)
			state_list.push_back(int8_t(*it));
		}

	// The compressIDLMatrix function first erases, then builds, both pattern_map and 
	// pattern_counts using the uncompressed data contained in mat
	std::vector<const CIPR_StateSet_t *> vecRowPtr;
	for (std::vector<CIPRStateVect>::const_iterator rowIt = m.matrix.begin(); rowIt != m.matrix.end(); ++rowIt)
		{
		const CIPRStateVect & v = *rowIt;
		vecRowPtr.push_back(&v[0]);
		}
	const CIPR_StateSet_t * const * matPtr = &vecRowPtr[0];
	num_patterns = compressIDLMatrix(m.ntaxa, m.nchar, matPtr);

	// size of likelihood_rate_site vector needs to be revisited if the number of rates subsequently changes 
	recalcRelativeRates();
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Builds `pattern_map' and `pattern_counts' using uncompressed data stored in `mat'.
*/
void TreeLikelihood::copyDataFromDiscreteMatrix(
  const CipresNative::DiscreteMatrix & mat)		/**< is the data source */
	{
	nTaxa = mat.getNTax();

	// The compressDataMatrix function first erases, then builds, both pattern_map and 
	// pattern_counts using the uncompressed data contained in mat
#if defined(INTERFACE_WITH_CIPRES)
	// changed signature because both compressDataMatrix and compressIDLMatrix now use 
	// buildPatternMapFromRawMatrix, which takes ntax and nchar as first two arguments
	num_patterns = compressDataMatrix(mat.getNTax(), mat.getNChar(), mat);
#else
	num_patterns = compressDataMatrix(mat);
#endif

	state_list = mat.getStateList(); 
	state_list_pos = mat.getStateListPos();

	// size of likelihood_rate_site vector needs to be revisited if the number of rates subsequently changes 
	recalcRelativeRates();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Builds `pattern_map' and `pattern_counts' using data stored in `sim_data'.
*/
void TreeLikelihood::copyDataFromSimData(
  SimDataShPtr sim_data)	/**< is the data source */
	{
	// Copy simulated data to pattern_map
	pattern_map = sim_data->getSimPatternMap();

	// Build up counts vector
	pattern_counts.clear();
	for (PatternMapType::iterator it = pattern_map.begin(); it != pattern_map.end(); ++it)
		{
		pattern_counts.push_back(it->second);
		}

	nTaxa = sim_data->getPatternLength();
	num_patterns = (unsigned)pattern_map.size();

	model->buildStateList(state_list, state_list_pos);

	// size of likelihood_rate_site vector needs to be revisited if the number of rates subsequently changes 
	recalcRelativeRates();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds data currently stored in `pattern_map' to the patterns already in `other'. Assumes that `pattern_length' for 
|	this SimData object is identical to the `pattern_length' of `other'.
*/
void TreeLikelihood::addDataTo(SimData & other)
	{
	if (pattern_map.empty())
		return;
	if (other.getTotalCount() == 0)
		{
		PHYCAS_ASSERT(nTaxa > 0);
		other.resetPatternLength(nTaxa);
		}
	PHYCAS_ASSERT(nTaxa == other.getPatternLength());
	for (PatternMapType::iterator it = pattern_map.begin(); it != pattern_map.end(); ++it)
		{
		PatternCountType count = it->second;
		VecStateList & other_pattern = other.getCurrPattern();
		std::copy(it->first.begin(), it->first.end(), other_pattern.begin());
		other.insertPattern(count);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls allocateInternalData() to add an InternalData structure to `nd' containing the conditional likelihood arrays
|	needed for likelihood calculations.
*/
void TreeLikelihood::prepareInternalNodeForLikelihood(
  TreeNode * nd)	/**< is the node to decorate */
	{
	//InternalData * ndID = nd->GetInternalData();
	if (!nd)
		{
		TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;
		InternalData * cl = allocateInternalData();
		nd->SetInternalData(cl, cl_deleter);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Decorates a tip node with data and add to the tree's StoreInternalNode vector.
*/
void TreeLikelihood::addDecoratedInternalNode(
  TreeShPtr t)		/**< is the tree to decorate */
	{
	TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;
	InternalData * cl = allocateInternalData();
    TreeNode * nd = t->AllocNewNode();
	nd->SetInternalData(cl, cl_deleter);
    t->StoreInternalNode(nd);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Decorates a tip node with data and add to the tree's tipStorage vector.
*/
void TreeLikelihood::addOrphanTip(
  TreeShPtr t,		/**< is the tree to decorate */
  unsigned row)		/**< is the row in the data matrix to associate with the added tip */
	{
	TreeNode::TipDataDeleter td_deleter	= &deallocateTipData;
	TipData * td = allocateTipData(row);
    TreeNode * nd = t->AllocNewNode();
	nd->SetTipData(td, td_deleter);
    nd->SetNodeNum(row);
    t->StoreLeafNode(nd);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Calls allocateTipData() for all tip nodes and allocateInternalData() for all internal nodes in the supplied Tree. 
|	Assumes each tip node number in the tree equals the appropriate row in the data matrix. 
*/
void TreeLikelihood::prepareForLikelihood( //POLBM TreeLikelihood::prepareForLikelihood
  TreeShPtr t)									/**< is the tree to decorate */
	{
	// If no_data is true, it means that calcLnL will always return 0.0 immediately and 
	// will thus never need the TipData or InternalData data structures
	if (no_data)
		return;

	TreeNode::TipDataDeleter		td_deleter	= &deallocateTipData;
	TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;

	// Put all existing conditional likelihood arrays already back into storage
	//storeAllCLAs(t);
	//cla_pool.clearStack();	//@POL this here only because prepareForLikelihood called in NCatMove::proposeNewState when ncat is increased

	for (preorder_iterator nd = t->begin(); nd != t->end(); ++nd)
		{
		if (nd->IsTip())
			{
			unsigned row = nd->GetNodeNumber();
			TipData * td = allocateTipData(row);
			nd->SetTipData(td, td_deleter);
			}
		else
			{
			InternalData * cl = allocateInternalData();
			nd->SetInternalData(cl, cl_deleter);
			}
		}

    //@POL  should decorate any nodes in tipStorage and nodeStorage now, but should wait until those 
    // are changed to vectors
	}

#if defined(INTERFACE_WITH_CIPRES)
/*----------------------------------------------------------------------------------------------------------------------
|	Helper function used by both compressDataMatrix and compressIDLMatrix that (re)builds `pattern_map', 
|	`pattern_counts' and `charIndexToPatternIndex' from a two-dimensional array `mat' of state codes (rows are taxa, 
|	columns are sites).
*/
unsigned TreeLikelihood::buildPatternMapFromRawMatrix(unsigned ntax, unsigned nchar, const CIPR_StateSet_t * const * mat)
	{
	pattern_map.clear();
	charIndexToPatternIndex.clear();

	typedef std::list<unsigned> IndexList;
	typedef std::map<VecStateList, IndexList> PatternToIndex;
	PatternToIndex patternToIndex;
	
	// Loop across each site in mat
	for (unsigned j = 0; j < nchar; ++j)
		{
		// Build up a vector representing the pattern of state codes at this site
		std::vector<int8_t> pattern;
		for (unsigned i = 0; i < ntax; ++i)
			{
			const CIPR_StateSet_t * row  = mat[i];
			const CIPR_StateSet_t   code = row[j];
			pattern.push_back(code);
			}

		// Add the pattern to the map if it has not yet been seen, otherwise increment 
		// the count of this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
		PatternMapType::iterator lowb = pattern_map.lower_bound(pattern);
		if (lowb != pattern_map.end() && !(pattern_map.key_comp()(pattern, lowb->first)))
			{
			// pattern is already in pattern_map, increment count
			lowb->second += 1.0;
			}
		else
			{
			// pattern has not yet been stored in pattern_map
			pattern_map.insert(lowb, PatternMapType::value_type(pattern, 1));
			}
		
		// Add the pattern to the map if it has not yet been seen, otherwise increment 
		// the count of this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
		PatternToIndex::iterator pToILowB = patternToIndex.lower_bound(pattern);
		if (pToILowB != patternToIndex.end() && !(patternToIndex.key_comp()(pattern, pToILowB->first)))
			pToILowB->second.push_back(j);
		else
			{
			IndexList ilist(1, j);
			patternToIndex.insert(pToILowB, PatternToIndex::value_type(pattern, ilist));
			}
		}

	// Copy counts to pattern_counts before returning
	pattern_counts.clear();
	pattern_counts.reserve(pattern_map.size());
	unsigned patternIndex = 0;
	for (PatternMapType::iterator mapit = pattern_map.begin(); mapit != pattern_map.end(); ++mapit, ++patternIndex)
		{
		pattern_counts.push_back(mapit->second);
		const IndexList & inds = patternToIndex[mapit->first];
		for (IndexList::const_iterator indIt = inds.begin(); indIt != inds.end(); ++indIt)
			charIndexToPatternIndex[*indIt] = patternIndex;
		}
	return (unsigned)pattern_map.size();
	}
#endif

#if defined(INTERFACE_WITH_CIPRES)
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
unsigned TreeLikelihood::compressIDLMatrix(unsigned ntax, unsigned nchar, const CIPR_StateSet_t * const * matrix)
	{
	return buildPatternMapFromRawMatrix(ntax, nchar, matrix);
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Copies data from `mat' to the map `pattern_map'. The resulting map holds pairs whose key is the pattern for one site
|	and whose value is a count of the number of sites having that pattern. The counts from `pattern_map' are transferred  
|	to the `pattern_counts' vector (vectors are more efficient containers for use during likelihood calculations).
*/
#if defined(INTERFACE_WITH_CIPRES)
unsigned TreeLikelihood::compressDataMatrix(unsigned ntax, unsigned nchar, const CipresNative::DiscreteMatrix & matrix) //POLBM TreeLikelihood::compressDataMatrix
	{
	const CIPR_StateSet_t * const * mat = matrix.getMatrix();
	return buildPatternMapFromRawMatrix(ntax, nchar, mat);
	}
#else
unsigned TreeLikelihood::compressDataMatrix(const CipresNative::DiscreteMatrix & mat) //POLBM TreeLikelihood::compressDataMatrix
	{
	pattern_map.clear();
	charIndexToPatternIndex.clear();
	unsigned ntax = mat.getNTax();
	unsigned nchar = mat.getNChar();

	typedef std::list<unsigned> IndexList;
	typedef std::map<VecStateList, IndexList> PatternToIndex;
	PatternToIndex patternToIndex;

	if (model->isCodonModel()) //@POL could move this test inside the taxon-loop to avoid being so redundant
		{
		// Loop across each triplet of sites in mat
		for (unsigned j = 0; j < nchar; j += 3)
			{
			// Suppose nchar=5 and j=3: not enough sites to make a second codon. In this example,
			// j+2 = 5, which equals nchar, so break out of loop over sites
			if (j+2 >= nchar)
				break;

			// Build up a vector representing the pattern of state codes at this site
			std::vector<int8_t> pattern;
			for (unsigned i = 0; i < ntax; ++i)
				{
				const int8_t * row  = mat.getRow(i);
				const int8_t   code1 = row[j];
				const int8_t   code2 = row[j+1];
				const int8_t   code3 = row[j+2];
				bool code1_ok = (code1 >= 0 && code1 < 4);
				bool code2_ok = (code2 >= 0 && code2 < 4);
				bool code3_ok = (code3 >= 0 && code3 < 4);
				if (code1_ok && code2_ok && code3_ok)
					{
                    const int8_t code = codon_state_codes[16*code1 + 4*code2 + code3]; // CGT = 27 = 16*1 + 4*2 + 3
					if (code > 60)
						throw XLikelihood(str(boost::format("Stop codon encountered for taxon %d at sites %d-%d") % (i+1) % (j+1) % (j+4)));
					pattern.push_back(code);
					}
				else
					pattern.push_back((int8_t)61);
				}

			//@POL below here same as the else block

			// Add the pattern to the map if it has not yet been seen, otherwise increment 
			// the count of this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
			PatternMapType::iterator lowb = pattern_map.lower_bound(pattern);
			if (lowb != pattern_map.end() && !(pattern_map.key_comp()(pattern, lowb->first)))
				{
				// pattern is already in pattern_map, increment count
				lowb->second += 1.0;
				}
			else
				{
				// pattern has not yet been stored in pattern_map
				pattern_map.insert(lowb, PatternMapType::value_type(pattern, 1));
				}
			
			// Add the pattern to the map if it has not yet been seen, otherwise increment 
			// the count of this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
			PatternToIndex::iterator pToILowB = patternToIndex.lower_bound(pattern);
			if (pToILowB != patternToIndex.end() && !(patternToIndex.key_comp()(pattern, pToILowB->first)))
				pToILowB->second.push_back(j);
			else
				{
				IndexList ilist(1, j);
				patternToIndex.insert(pToILowB, PatternToIndex::value_type(pattern, ilist));
				}
			}
		}
	else
		{
		// Loop across each site in mat
		for (unsigned j = 0; j < nchar; ++j)
			{
			// Build up a vector representing the pattern of state codes at this site
			std::vector<int8_t> pattern;
			for (unsigned i = 0; i < ntax; ++i)
				{
				const int8_t * row  = mat.getRow(i);
				const int8_t   code = row[j];
				pattern.push_back(code);
				}

			// Add the pattern to the map if it has not yet been seen, otherwise increment 
			// the count of this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
			PatternMapType::iterator lowb = pattern_map.lower_bound(pattern);
			if (lowb != pattern_map.end() && !(pattern_map.key_comp()(pattern, lowb->first)))
				{
				// pattern is already in pattern_map, increment count
				lowb->second += 1.0;
				}
			else
				{
				// pattern has not yet been stored in pattern_map
				pattern_map.insert(lowb, PatternMapType::value_type(pattern, 1));
				}
			
			// Add the pattern to the map if it has not yet been seen, otherwise increment 
			// the count of this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
			PatternToIndex::iterator pToILowB = patternToIndex.lower_bound(pattern);
			if (pToILowB != patternToIndex.end() && !(patternToIndex.key_comp()(pattern, pToILowB->first)))
				pToILowB->second.push_back(j);
			else
				{
				IndexList ilist(1, j);
				patternToIndex.insert(pToILowB, PatternToIndex::value_type(pattern, ilist));
				}
			}
		}

	// Copy counts to pattern_counts before returning
	pattern_counts.clear();
	pattern_counts.reserve(pattern_map.size());
	unsigned patternIndex = 0;
	for (PatternMapType::iterator mapit = pattern_map.begin(); mapit != pattern_map.end(); ++mapit, ++patternIndex)
		{
		pattern_counts.push_back(mapit->second);
		const IndexList & inds = patternToIndex[mapit->first];
		for (IndexList::const_iterator indIt = inds.begin(); indIt != inds.end(); ++indIt)
			charIndexToPatternIndex[*indIt] = patternIndex;
		}

	return (unsigned)pattern_map.size();
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a string representation of the supplied `state'. For example, if `state' equals 1 (where model is a 
|	standard DNA model), the string returned would be "C". If, however, `state' was 7 (again, standard DNA model), then
|	there is some ambiguity present and the string returned might look something like "{AC}" (the string returned 
|	depends of course on the actual meaning of the global state code 7).
*/
std::string TreeLikelihood::getStateStr(
  int8_t state)	const /**< is the global state code to be converted to a std::string */
	{
	std::string s;
	const int8_t nsPlusOne = num_states + 1;

	if (state < nsPlusOne)
		{
		// either no ambiguity or complete ambiguity
		s << model->lookupStateRepr((int)state);
		}
	else
		{
		// `state' represents partial ambiguity

		// First, find location of the definition of `state' in the global state list
		unsigned pos = state_list_pos[(unsigned)state];
		VecStateList::const_iterator it = state_list.begin() + pos;

		// Now get the number of basic states composing `state'
		unsigned n = *it++;

		// Walk down global state list converting states into strings
		s << "{";
		for (unsigned i = 0; i < n; ++i)
			{
			s << model->lookupStateRepr((int)*it++);
			}
		s << "}";
		}
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets private data member `nevals' to 0.
*/
void TreeLikelihood::resetNEvals()
	{
	nevals = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor providing access to the current value of the private data member `nevals'.
*/
unsigned TreeLikelihood::getNEvals()
	{
	return nevals;
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Sets the boolean data member `debugging_now' to the value specified.
*/
void TreeLikelihood::setDebug(
  bool on)  /**< is the value to which `debugging_now' should be set */
    {
    debugging_now = on;
    }
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Assuming TreeLikelihood::compressDataMatrix has been called, so that `pattern_map' is up-to-date, returns a string 
|	listing all observed patterns and their frequencies.
*/
std::string TreeLikelihood::listPatterns(
  bool show_coded_states)	/**< if true, output global state codes used internally; otherwise, do the translation back to the original state representations */
	{
	std::vector<std::string> state_repr;
	std::string s;
	unsigned i = 0;
	PatternMapType::iterator it = pattern_map.begin();
	for (; it != pattern_map.end(); ++it, ++i)
		{
		const std::vector<int8_t> & p = it->first;
		PatternCountType c = it->second;
		s << str(boost::format("%6d %6.1f ") % i % c);
		unsigned ntax = (unsigned)p.size();
		if (show_coded_states)
			{
			for (unsigned j = 0; j < ntax; ++j)
				{
				s << str(boost::format("%d ") % (int)(p[j]));
				}
			}
		else
			{
			for (unsigned j = 0; j < ntax; ++j)
				{
				s << str(boost::format("%s") % getStateStr(p[j]));
				}
			}
		s << '\n';
		}
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First erases `pattern_repr', then builds it into a vector of pattern representation strings. When the function 
|	returns, each element of `pattern_repr' is a string containing global state codes separated by whitespace. This 
|	function is intended to be used for debugging purposes, where it is sometimes helpful to be sure of exactly which
|	data pattern is being operated upon. Thus, no particular effort has been made to make this function efficient.
*/
void TreeLikelihood::buildPatternReprVector(std::vector<std::string> & pattern_repr, TreeShPtr t)
	{
	unsigned nTips		= t->GetNTips();
	//unsigned nStates	= getNStates();
	//unsigned nPatterns	= getNPatterns();

	pattern_repr.clear();
		pattern_repr.reserve(num_patterns);

	std::cerr << "\nglobal_pos:" << std::endl;
	for (unsigned z = 0; z < state_list_pos.size(); ++z)
		{
		std::cerr << str(boost::format("%3d") % z) << "  " << (int)state_list_pos[z] << std::endl;
		}

	// Recreate the data matrix in terms of tip-specific state codes
	int8_t * * m = NewTwoDArray<int8_t>(nTips, num_patterns);

	for (preorder_iterator nd = t->begin(); nd != t->end(); ++nd)
		{
		if (nd->IsTip())
			{
			// get node number i
			unsigned i = nd->GetNodeNumber();
			PHYCAS_ASSERT(i >= 0 && i < nTips);

			// get tip-specific state code array
			TipData * tipData = nd->GetTipData();
			PHYCAS_ASSERT(tipData);
			const int8_t * tipCodes = tipData->getConstStateCodes();

			// get position vector that allows us to translate tip-specific state codes into global state codes
			const std::vector<unsigned int> & global_statelist_pos = tipData->getConstStateListPos();

			std::cerr << "\nglobal_statelist_pos vector for node " << i << ": " << std::endl;
			for (unsigned z = 0; z < global_statelist_pos.size(); ++z)
				{
				std::cerr << "  " << (int)global_statelist_pos[z] << std::endl;
				}

			// fill row i of data matrix
			for (unsigned j = 0; j < num_patterns; ++j)
				{
				int8_t global_code = tipCodes[j];
				int8_t offset = global_code - ((int8_t)num_states + 1);
				if (offset >= 0)
					{
					unsigned pos = global_statelist_pos[offset];
					for (unsigned m = 0; m < (unsigned)state_list_pos.size(); ++m)
						{
						if (state_list_pos[m] == pos)
							global_code = (int8_t)(m);
						}
					}
				m[i][j] = global_code;

				//std::cerr << j << " | " << (int)tipCodes[j] << " | ";
				//if (offset >= 0)
				//	std::cerr << (int)global_statelist_pos[offset];
				//else
				//	std::cerr << "(" << (int)offset << ")";
				//std::cerr << std::endl;
				}
			}
		}

	for (unsigned j = 0; j < num_patterns; ++j)
		{
		std::string s;
		for (unsigned i = 0; i < nTips; ++i)
			{
			s << (int)m[i][j] << " ";
			}
		pattern_repr.push_back(s);
		}

	DeleteTwoDArray<int8_t>(m);
	}

}	// namespace phycas
