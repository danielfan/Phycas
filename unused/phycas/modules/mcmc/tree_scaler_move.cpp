#include "phycas/force_include.h"
#include "phycas/rand/lot.hpp"
#include "phycas/trees/tree.hpp"
#include "phycas/trees/tree_node.hpp"
#include "phycas/modules/mcmc/tree_scaler_move.hpp"
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::log;
	using std::exp;
#endif


/*--------------------------------------------------------------------------------------------------------------------------
|	Constructor passes `t' to the base class TreeMove constructor, sets `lambda' and `scalingFactor' both to 1.0.
*/
TreeScalerMove::TreeScalerMove(
  Tree *t)	/**< is the pointer to the tree, which is passed along to base class TreeMove constructor */
	: TreeMove(t)
	{
	assert(t != NULL);
	scalingFactor = 1.0;
	lambda = 1.0;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Nothing to be done.
*/
TreeScalerMove::~TreeScalerMove()
	{
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is accepted.
*/
void TreeScalerMove::Accept()
	{
	scalingFactor = 1.0;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Causes tree to be returned to its state just prior to proposing the move.
*/
void TreeScalerMove::Revert()
	{
	for (TreeNode *nd = tree->GetFirstPostorder(); nd != NULL; nd = nd->GetNextPostorder())
		{
		double originalEdgeLen = nd->GetFltEdgeLen()/scalingFactor;
		nd->SetFltEdgeLen(originalEdgeLen);
		nd->InvalidateAttrDown(true, nd->GetParent());
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move. The Hastings ratio is scalingFactor^{n}, where n is the 
|	number of edges in the tree. Assumes tree is unrooted.
*/
double TreeScalerMove::GetLnHastingsRatio() const
	{
	// n is the number of edges in the tree
	double n = (double)(tree->GetNNodes() - 1);	//@ assumes tree is unrooted at the moment
	return n*log(scalingFactor);	
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	This move does not change the model dimension, so the Jacobian is irrelevant.
*/
double TreeScalerMove::GetLnJacobian() const
	{
	return 0.0;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Multiply the lengths of all edges in the tree by a common `scalingFactor', which is determined as follows:
|>
|	`scalingFactor' = exp(lambda*(r.Uniform() - 0.5))
|<
*/
void TreeScalerMove::ProposeNewState()
	{
	// Make sure MCMCMove::SetRandomNumberGenerator has been called to set rng
	//
	assert(rng != NULL); 

	// Select a new scalingFactor
	//
	scalingFactor = exp(lambda*(rng->Uniform() - 0.5));

	for (TreeNode *nd = tree->GetFirstPostorder(); nd != NULL; nd = nd->GetNextPostorder())
		{
		double newEdgeLen = nd->GetFltEdgeLen()*scalingFactor;
		nd->SetFltEdgeLen(newEdgeLen);
		nd->InvalidateAttrDown(true, nd->GetParent());
		}
	}
