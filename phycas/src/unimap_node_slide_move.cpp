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
#include <fstream>
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/basic_tree_node.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/unimap_node_slide_move.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/internal_data.hpp"
#include "phycas/src/tip_data.hpp"

#include "boost/format.hpp"

namespace phycas
{



// Choose random internal node origNode and randomly choose one of origNode's
//	children to call x (the other is y)
// Here is the actual layout of the tree in memory
//
//      x  y        x = child that is swapped
//       \/         y = x's sibling
//   z   origNode	z = origNode's sibling
//    \ /           w = origNodePar's parent
//    origNodePar
//    /
//   w
// the likelihood calculations are done assuming that the tree is ((a,b),(c,d)) with the internal length from origNd)
void UnimapNodeSlideMove::ProposeStateWithTemporaries(ChainManagerShPtr & p)
	{
	PHYCAS_ASSERT(false);
	const double u = rng->Uniform(FILE_AND_LINE);
	if (u < 0.5)
		{
		if (u < 0.25)
			{
			moveType = MOVE_X_TOWARD_W;
			nniPartOfPath = prev_ndP_len;
			}
		else
			{
			moveType = moveType;
			nniPartOfPath = prev_z_len;
			}
		startingDist = nniPartOfPath + origNode->GetEdgeLen();
		}
	else
		{
		if (u < 0.75)
			{
			moveType = MOVE_W_TOWARD_Y;
			startingDist = prev_z_len;
			}
		else
			{
			moveType = MOVE_Z_TOWARD_Y;
			startingDist = prev_ndP_len;
			}
		nniPartOfPath = prev_nd_len + startingDist;
		}
	const double pathToSlideAlong = prev_y_len + startingDist;


	PHYCAS_ASSERT(a == y);
	PHYCAS_ASSERT(b == x);
	PHYCAS_ASSERT(c == w);
	PHYCAS_ASSERT(d == z);
	PHYCAS_ASSERT(aLenNd == a);
	PHYCAS_ASSERT(bLenNd == b);
	PHYCAS_ASSERT(cLenNd == origNodePar);
	PHYCAS_ASSERT(dLenNd == d);

	// choose a distance to displace the internal node, add it to the path, and then reflect to keep it legal...
	double displacementWindow = windowSize*(rng->Uniform(FILE_AND_LINE) - 0.5);
	

	finalDist = startingDist + displacementWindow;
	while (finalDist < 0.0 || finalDist > pathToSlideAlong)
		{
		if (finalDist < 0.0)
			finalDist = -finalDist;
		if (finalDist > pathToSlideAlong)
			finalDist = (2.0*pathToSlideAlong) - finalDist;
		}
	
	
	if (moveType == MOVE_X_TOWARD_W || moveType == MOVE_X_TOWARD_Z)
		{
		if (finalDist < nniPartOfPath)
			{
			if (moveType == MOVE_X_TOWARD_W)
				{ // z and y will be sisters
				// this makes b = z,   c = w,  and   d = x;
				std::swap(b, d);
				std::swap(bLenNd, dLenNd);
				std::swap(bTipData, dTipData);
				cLenNd->SetEdgeLen(finalDist);
				}
			else
				{ // w and y will be sisters
				// this makes b = w,   c = x,  and   d = z;
				std::swap(b, c);
				std::swap(bLenNd, cLenNd);
				std::swap(bTipData, cTipData);
				dLenNd->SetEdgeLen(finalDist);
				}
			origNode->SetEdgeLen(nniPartOfPath - finalDist);
			a->SetEdgeLen(pathToSlideAlong - nniPartOfPath);
			}
		else
			{ // no topology change
			a->SetEdgeLen(pathToSlideAlong - finalDist);
			origNode->SetEdgeLen(finalDist - nniPartOfPath);
			}
		}
	else
		{
		if (finalDist > nniPartOfPath)
			{
			if (moveType == MOVE_Z_TOWARD_Y)
				{ // z and y will be sisters
				// this makes b = z,   c = w,  and   d = x;
				std::swap(b, d);
				std::swap(bLenNd, dLenNd);
				std::swap(bTipData, dTipData);
				cLenNd->SetEdgeLen(nniPartOfPath);
				}
			else
				{ // w and y will be sisters
				// this makes b = w,   c = x,  and   d = z;
				std::swap(b, c);
				std::swap(bLenNd, cLenNd);
				std::swap(bTipData, cTipData);
				dLenNd->SetEdgeLen(nniPartOfPath);
				}
			origNode->SetEdgeLen(finalDist - nniPartOfPath);
			a->SetEdgeLen(pathToSlideAlong - finalDist);
			}
		else
			{ // no topology change
			if (moveType == MOVE_Z_TOWARD_Y)
				cLenNd->SetEdgeLen(finalDist);
			else
				dLenNd->SetEdgeLen(finalDist);
			origNode->SetEdgeLen(nniPartOfPath - finalDist);
			}
		}
	
	PHYCAS_ASSERT(aLenNd == a);
	PHYCAS_ASSERT((bLenNd == b && cLenNd == origNodePar) || (bLenNd == origNodePar && cLenNd == c));
	PHYCAS_ASSERT(dLenNd == d);
        
	curr_ln_prior = calcEdgeLenLnPrior(*a, aLenNd->GetEdgeLen(), p)
					+ calcEdgeLenLnPrior(*b, bLenNd->GetEdgeLen(), p)
					+ calcEdgeLenLnPrior(*c, cLenNd->GetEdgeLen(), p)
					+ calcEdgeLenLnPrior(*d, dLenNd->GetEdgeLen(), p)
					+ calcEdgeLenLnPrior(*origNode, origNode->GetEdgeLen(), p);
	}


/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void UnimapNodeSlideMove::accept()
	{
	++nMovesAccepted;

	PHYCAS_ASSERT(origNode);
	PHYCAS_ASSERT(origNode->GetParent() == origNodePar);

	if (moveType == MOVE_X_TOWARD_W || moveType == MOVE_X_TOWARD_Z)
		{
		if (finalDist < nniPartOfPath)
			{
			if (moveType == MOVE_X_TOWARD_W)
				tree_manipulator.NNISwap(z, x);
			else
				tree_manipulator.NNISwap(z, y);
			}
		}
	else
		{
		if (finalDist > nniPartOfPath)
			{
			if (moveType == MOVE_Z_TOWARD_Y)
				tree_manipulator.NNISwap(z, x);
			else
				tree_manipulator.NNISwap(z, y);
			}
		}
	
	if (doSampleInternalStates)
		resampleInternalNodeStates(post_root_posterior->getCLA(), post_cla->getCLA());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The default constructor sets `lambda' to the default value (0.2), sets `topol_changed' to false, and `m' and `mstar'
|	to 0.0. All other data members are automatically initialized (shared pointers) or are initialized via a call to 
|	reset().
*/
UnimapNodeSlideMove::UnimapNodeSlideMove() : UnimapTopoMove(),
  moveType(MOVE_X_TOWARD_W),
  windowSize(0.05)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move.
*/
double UnimapNodeSlideMove::getLnHastingsRatio() const
	{
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Jacobian for this move.
*/
double UnimapNodeSlideMove::getLnJacobian() const
	{
	return 0.0;
	}

}	// namespace phycas
