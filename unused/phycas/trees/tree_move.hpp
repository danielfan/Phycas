#ifndef PHO_TREEMOVE_H
#define PHO_TREEMOVE_H

#include "phycas/rand/lot.hpp"
#include "phycas/modules/mcmc/mcmc_move.hpp"
#include "phycas/trees/tree_manip.hpp"

#include "phycas/trees/tree.hpp"	//inlines only

//	Base class for mcmc tree moves
class TreeMove : public MCMCMove , public TreeManip
	{
	protected :
		TreeNode	*ndBase;	/* most ancestral node involved in the move, used as the center of the likelihood calcuation (and in Revert) */

	public	:
					TreeMove(Tree *);
					virtual ~TreeMove();

		TreeNode	*GetAnAffectedNode();
		bool		IsMetropolisHastings() const;
	};
	
inline TreeMove::TreeMove(
  Tree *t)
  : TreeManip(t)
	{
	assert(t != NULL);
	ndBase = t->GetRoot();
	}

inline TreeMove::~TreeMove()
	{
	}

inline TreeNode *TreeMove::GetAnAffectedNode()
	{
	return ndBase;
	}

inline bool TreeMove::IsMetropolisHastings() const
	{
	return true;
	}

#endif

