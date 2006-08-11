//#define POL_MESQUITE_COLORING

#include "phycas/force_include.h"
#include "phycas/modules/mcmc/larget_simon_move.hpp"
#include "phycas/trees/tree_node.hpp"
#include "phycas/trees/tree.hpp"
#include "phycas/rand/lot.hpp"
typedef unsigned nTaxaInt;

#define KEEP_BRANCHES_LEGAL
#define MAX_LEGAL_BRLEN 285

#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::log;
	using std::pow;
	using std::exp;
#endif

/*--------------------------------------------------------------------------------------------------------------------------
|	The default constructor sets lambda to the supplied value `lam' and initializes all other data members to NULL (if 
|	pointer) or 0.0 (if double valued). The `topologyChanged' data member is initialized to false.
*/
LargetSimonMove::LargetSimonMove(
  Tree *t,		/* passed along to base class TreeMove constructor */
  double lam)	/* value of tuning parameter lambda */
	: TreeMove(t)
	{
	topologyChanged	= false;
	lambda			= lam;
	m				= 0.0;
	mstar			= 0.0;
	Reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Performs a local perturbation using a generalization of the algorithm "LOCAL Without a Molecular Clock" described by 
|	Larget and Simon (1999. Mol. Biol. Evol. 16(6): 750-759). This version allows polytomies (except the case of the star
|	tree, for which the required three-continuous-edge segment cannot be identified).
|	
|	     x  y  z
|	      \ | /
|          \|/
|	  b  c  v
|	   \ | /
|	    \|/
|	     u
|	     |
|	     a <-- may or may not be leaf node at which tree is rooted
|	
|	Pick a random interior node v whose parent is not the root (i.e. avoid the internal node directly connected to the leaf
|	node at which the tree is rooted. Let u be the parent of v. Let a be a randomly-chosen child of u (note that in this case
|	u's parent is considered a "child" of u. In the figure above, we (by chance) chose the parent of u to be a, but we could
|	have chosen any of u's real children (except v) to be the "a" node.. Let b, c, ..., be the other "children" of u, not 
|	including v. Let x be a randomly chosen child of v (and here a child is really a child). Let y, z, ..., be the other
|	children of v.
|	
|	      b   c   y   z
|	       \ /     \ /
|	a ===== u ===== v ===== x
|	
|	The path represented by the double line above is either shrunk or expanded by a factor m*, where
|	
|	m* = m*exp(lambda*(r.Uniform() - 0.5))
|	
|	Then, one of {u, v} is chosen at random to move. Let's say for illustration that u was chosen. u is moved (along with b
|	and c) to a random point along the main path from a to x. If this makes the node across v, then the equivalent of an NNI 
|	rearrangement will be effected (a swapped with x); otherwise, the move will only require adjusting branch lengths.
*/
void LargetSimonMove::ProposeNewState()
	{
	double x, y, z, xstar;

	// Make sure MCMCMove::SetRandomNumberGenerator has been called to set rng
	//
	assert(rng != NULL); 

	// Begin by resetting all the data members involved with reverting a move
	//
	Reset();
	topologyChanged	= false;

	// Must avoid first internal node (only child of root), so number of acceptable nodes
	// is one fewer than the number of internal nodes
	//
	nTaxaInt numAcceptableNodes = tree->GetNInternals() - 1;

	// Select an internal node whose parent is not the root node to serve as ndY, whose branch will form the middle
	// segment in the path of three contiguous segments to be modified.
	//
	nTaxaInt ypos = rng->SampleUInt(numAcceptableNodes);
	nTaxaInt i = 0;
	ndY = NULL;
	for (ndY = tree->GetFirstPreorder(); ndY != NULL; ndY = ndY->GetNextPreorder())
		{
		if (!ndY->IsLeafOrRoot() && !ndY->GetParent()->IsRoot())
			{
			if (i == ypos)
				break;
			++i;
			}
		}
	assert(ndY->GetLeftChild() != NULL);
	assert(ndY->GetParent() != NULL);
	assert(!ndY->GetParent()->IsRoot());

	// Save ndY's edge length in case revert is needed.
	//
	origY = ndY->GetFltEdgeLen();

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
	origX = ndX->GetFltEdgeLen();

	// Set ndZ equal to a randomly-chosen child of U = ndY->par, ignoring ndY but treating U->par as if it were a child
	// If U->par is chosen to be ndZ, let ndZ equal U instead because the three nodes (ndX, ndY and ndZ) should be the nodes
	// that manage the relevant edge lengths in the move.
	//
	TreeNode *U = ndY->GetParent();
	assert(U != NULL);
	unsigned uchildren = U->CountChildren();
	which_child = rng->SampleUInt(uchildren);
	if (which_child == 0)
		{
		// Selected "child" is actually U's parent
		//
		ndZ = U;
		origZ = U->GetFltEdgeLen();

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
		origZ = ndZ->GetFltEdgeLen();

		//	Set ndBase to the deepest affected node
		//
		ndBase = U;
		}

	m = origX + origY + origZ;
	mstar = m*exp(lambda*(rng->Uniform() - 0.5));
	x = origX*mstar/m;
	y = origY*mstar/m;
	z = origZ*mstar/m;

	xstar = rng->Uniform() * mstar;

	// Decide whether to move ndY down or node U up
	//
	bool moving_Y = true;
	if (rng->Uniform() < 0.5)
		moving_Y = false;
	bool moving_U = !moving_Y;

#if 0 //TESTING(LargetSimonMove::ProposeNewState)
	if (moving_Y)
		xstar = x + y + z/2.0;
	else
		xstar = x/2.0 + y + z;
#endif

	if (moving_Y && (xstar <= x + y))
		{
		// No change in topology
		//
		ndX->SetFltEdgeLen(xstar);
		ndY->SetFltEdgeLen(x + y - xstar);
		ndZ->SetFltEdgeLen(z);
		}
	else if (moving_U && (xstar <= y + z))
		{
		// No change in topology
		//
		ndX->SetFltEdgeLen(x);
		ndY->SetFltEdgeLen(y + z - xstar);
		ndZ->SetFltEdgeLen(xstar);
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
			NNISwap(swap1,swap2);
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

			NNISwapSpecial(swap1);

#if 0	// this code now in NNISwapSpecial
			TreeNode *R = ndZ->FindRightChild();
			if (R == ndY)
				R = ndX;
			TreeNode *L = ndZ->GetLeftChild();

			while (ndZ->GetLeftChild() != ndY)
				LChildToLSib(ndZ, ndX);
			while (ndY->GetRightSib() != NULL)
				RChildToRSib(ndZ, ndX);
			if (ndY != L)
				{
				while (ndY->GetLeftChild() != L)
					LChildToLSib(ndY, ndY);
				}
			while (R->GetRightSib() != NULL)
				RChildToRSib(ndY, ndY);
#endif
			}

		if (moving_Y)
			{
			ndX->SetFltEdgeLen(x + y);
			ndY->SetFltEdgeLen(xstar - x - y);
			ndZ->SetFltEdgeLen(x + y + z - xstar);
			}
		else
			{
			ndX->SetFltEdgeLen(x + y + z - xstar);
			ndY->SetFltEdgeLen(xstar - y - z);
			ndZ->SetFltEdgeLen(y + z);
			}

		NXS_ASSERT(ndX->GetFltEdgeLen() > 0.0);
		NXS_ASSERT(ndY->GetFltEdgeLen() > 0.0);
		NXS_ASSERT(ndZ->GetFltEdgeLen() > 0.0);

		topologyChanged = true;
		}

#	if defined DOUBLE_ATTRS
		ndX->InvalidateAttrDown(true, ndBase);
		ndY->InvalidateAttrDown(true, ndBase);
		ndZ->InvalidateAttrDown(true, ndBase);
		/*
		if (swap1)
			swap1->InvalidateAttrDown(false,ndBase);
		if (swap2)
			swap2->InvalidateAttrDown(false,ndBase);
		*/
#	else
		ndX->InvalidateAttrDown(true, ndBase);
		ndY->InvalidateAttrDown(true, ndBase);
		ndZ->InvalidateAttrDown(true, ndBase);
		/*
		if (swap1)
			swap1->InvalidateAttrDown(false,ndBase);
		if (swap2)
			swap2->InvalidateAttrDown(false,ndBase);
		*/
#	endif

#	if defined(POL_MESQUITE_COLORING)
	tree->RefreshSelectionStatus();
	ndX->SetSelectionStatus(true);
	ndY->SetSelectionStatus(true);
	ndZ->SetSelectionStatus(true);
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reverses move made in ProposeNewState. Assumes ndX, ndY, and ndZ are non-NULL, which will be true if 
|	ProposeNewState was just called.
*/
void LargetSimonMove::Revert()
	{
	assert(ndX != NULL);
	assert(ndY != NULL);
	assert(ndZ != NULL);
	assert(topologyChanged ? (swap1 != NULL && swap2 != NULL) : (swap1 == NULL && swap2 == NULL));

	if (topologyChanged)
		{
		if (swap2 == ndZ)
			{
			// If swap2 equals ndZ, then swap2 was a child of ndBase and we were able to use
			// the standard NNISwap function to swap the two nodes
			//
			NNISwap(swap1,swap2);
			}
		else 
			{
			// If swap2 is ndZ's parent, then swap2 is ndBase (i.e. it is the "child" node below the
			// lower of the two adjacent internal nodes involved in the swap) and we had to use the
			// NNISwapSpecial function to perform the rearrangment
			//
			NNISwapSpecial(swap1);
			}
		}

	ndX->SetFltEdgeLen(origX);
	ndY->SetFltEdgeLen(origY);
	ndZ->SetFltEdgeLen(origZ);

#	if defined DOUBLE_ATTRS
//@POL Mark, not sure this is correct anymore. Let me explain the system I'm using now before trying to rewrite this section.
		ndX->InvalidateAttrDownRejecting(true, ndBase);
		ndY->InvalidateAttrDownRejecting(true, ndBase);
		ndZ->InvalidateAttrDownRejecting(true, ndBase);
	for (TreeNode *nd = tree->GetFirstPreorder(); nd; nd = nd->nextPreorder)
		{
		if (nd->IsTip() == false)
			{
			for(int j = 0; j < tree->nInternalAttr; ++j)
				{
				if (nd->otherAttr[j])
					{
					TreeNodeAttr* temp = nd->attr[j];
					nd->attr[j] = nd->otherAttr[j];
					nd->otherAttr[j] = temp;
					}
				}
			}
		}		
#	else
		ndX->InvalidateAttrDown(true,ndBase);
		ndY->InvalidateAttrDown(true,ndBase);
		ndZ->InvalidateAttrDown(true,ndBase);
#	endif
	
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted only does anything if multiple copies of the attributes are being stored (so old ones can
|	be marked as "dead", and the current Attr can be "rotated" (so that the last move is attributes are always in otherAttr 
|	not attr place
*/
void LargetSimonMove::Accept()
	{
#	if defined DOUBLE_ATTRS
		ndX->InvalidateAttrRecursivelyAccepting(true);
		ndY->InvalidateAttrRecursivelyAccepting(true);
		ndZ->InvalidateAttrRecursivelyAccepting(true);
	for (TreeNode *nd = tree->GetFirstPreorder(); nd; nd = nd->nextPreorder)
		{
		if (nd->IsTip() == false)
			{
			for(int j = 0; j < tree->nInternalAttr; ++j)
				{
				if (nd->otherAttr[j])
					{
					TreeNodeAttr* temp = nd->attr[j];
					nd->attr[j] = nd->otherAttr[j];
					nd->otherAttr[j] = temp;
					}
				}
			}
		}		
#	endif
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Forgets information saved to enable reverting a proposed move. Like the Forget function, except the Forget function
|	also calls the base class Forget function and thus forgets more information than does Reset.
*/
inline void LargetSimonMove::Reset()
	{
	swap1		= NULL;
	swap2		= NULL;
	ndX			= NULL;
	ndY			= NULL;
	ndZ			= NULL;
	ndBase		= NULL;
	origX		= 0.0;
	origY		= 0.0;
	origZ		= 0.0;
	}



#if 0 //HIDE_PROTECTION_FUNCS
	LargetSimonMove &LargetSimonMove::operator=(const LargetSimonMove &)
		{
		assert(0);
		throw BadEqOpXcp("LargetSimonMove");
		}
#endif

