#include "phycas/force_include.h"
#include "phycas/modules/mcmc/edge_move.hpp"
#include "phycas/rand/lot.hpp"
#include "phycas/trees/tree.hpp"
#include "phycas/trees/tree_node.hpp"
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::log;
	using std::pow;
	using std::exp;
#endif
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
EdgeMove::EdgeMove(
  Tree *t)	/* passed along to base class TreeMove constructor */
	: TreeMove(t)
	{
	assert(t != NULL);
	origEdgelen.f = 0.0;
	origNode = NULL;
	lambda = 1.0;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Nothing to be done.
*/
EdgeMove::~EdgeMove()
	{
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is accepted.
*/
void EdgeMove::Accept()
	{
	origEdgelen.f = 0.0;
	origNode = NULL;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Causes tree to be returned to its state just prior to proposing the move.
*/
void EdgeMove::Revert()
	{
	origNode->SetFltEdgeLen(origEdgelen.f);
	origEdgelen.f = 0.0;
	origNode->InvalidateAttrDown(true, origNode->GetParent());
	origNode = NULL;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Inlined function that returns the natural log of the Hastings ratio for this move. The Hastings ratio is:
|>
|	Pr(new state -> old state)   [1/(lambda m)]
|	-------------------------- = --------------- = m* / m
|	Pr(old state -> new state)   [1/(lambda m*)]
|>
|	Here, m* is the new length of the modified edge and m is the length before the move is proposed.
*/
double EdgeMove::GetLnHastingsRatio() const
	{
	return log(origNode->GetFltEdgeLen()) - log(origEdgelen.f);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	This move does not change the model dimension, so the Jacobian is irrelevant.
*/
double EdgeMove::GetLnJacobian() const
	{
	return 0.0;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Chooses a random edge and changes its present length m to a new length m* using the formula, where lambda is a tuning
|	parameter.
|>
|	m* = m*exp(lambda*(r.Uniform() - 0.5))
|>
*/
void EdgeMove::ProposeNewState()
	{
	// Make sure MCMCMove::SetRandomNumberGenerator has been called to set rng
	//
	assert(rng != NULL); 

	// Choose edge randomly.
	//
	unsigned numEdges = tree->GetNNodes() - 1;
	unsigned k = rng->SampleUInt(numEdges);
	unsigned i = 0;
	for (origNode = tree->GetFirstPreorder(); origNode != NULL; origNode = origNode->GetNextPreorder())
		{
		// All nodes have an edge associated with them except for the root
		//
		if (!origNode->IsRoot())
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
	double m = origNode->GetFltEdgeLen();
	double mstar = m*exp(lambda*(rng->Uniform() - 0.5));
	origNode->SetFltEdgeLen(mstar);
	origNode->InvalidateAttrDown(true, origNode->GetParent());
	}

