#include "phycas/force_include.h"
#include "phycas/modules/mcmc/move_schedule.hpp"
#include "phycas/rand/lot.hpp"
#include "phycas/modules/mcmc/bush_master.hpp"
#include "phycas/modules/mcmc/edge_move.hpp"
#include "phycas/modules/mcmc/larget_simon_move.hpp"
#include "phycas/modules/mcmc/tree_scaler_move.hpp"
#include "phycas/command/mcmc_settings.hpp"


void MoveSchedule::ResetMoves(Tree * tree, LotShPtr r)
	{
	lsm	= LargetSimonMoveShPtr(new LargetSimonMove(tree, 0.2));
	bushMaster = BushMasterShPtr(new BushMaster(r, tree));
	edgeMove = EdgeMoveShPtr(new EdgeMove(tree));
	treeScalerMove = TreeScalerMoveShPtr(new TreeScalerMove(tree));
	
	lsm->SetRandomNumberGenerator(r);
    edgeMove->SetRandomNumberGenerator(r);
	treeScalerMove->SetRandomNumberGenerator(r);
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Builds `vProbTopoMove' vector from weights provided by user. Afterwards, vProbTopoMove[i] contains the probability
|	that the ith topology move will be attempted in any given MCMC step.
*/
void MoveSchedule::Reset(const MCMCSettings & settings)
	{
	doGibbs = settings.gibbsEvery > 0;
	gibbsEvery = settings.gibbsEvery;
	vProbTopoMoves.clear();
	const double sum = settings.localMoveWeight + settings.bushMoveWeight + settings.scalerMoveWeight;
		// It is possible for user to want to perform only Gibbs sampling of model parameters
		// on the starting tree, in which case the weights for MCMC moves will sum to zero and
		// vProbTopoMove will remain empty
	if (sum > 0.0)
		{
		vProbTopoMoves.push_back(settings.localMoveWeight/sum);
		vProbTopoMoves.push_back(settings.bushMoveWeight/sum);
		vProbTopoMoves.push_back(settings.scalerMoveWeight/sum);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns an index that indicates which tree move to attempt next. Walks through vProbTopoMove vector until cumulative
|	probability exceeds random uniform deviate drawn beforehand.
*/
MCMCMove * MoveSchedule::GetNextTreeMove(Lot & r, const Tree & tree) const
	{
	assert(!vProbTopoMoves.empty());
	//@POL, this code should be equivalent to the loop below, 
	//  but SingleMulitinomial returns an index in the range  [0, vProbTopoMove.size() ), 
	//  not [0, vProbTopoMove.size() ].  This is OK, right? (I'm assuming that return sz was unreachable,
	//  because vProbTopoMove sum to 1).
	const unsigned m = r.SingleMulitinomial(&vProbTopoMoves[0], (unsigned)vProbTopoMoves.size()); 
	if (m == MoveSchedule::TreeMove_Bush)
		return bushMaster.get();
	if (m == MoveSchedule::TreeMove_Local) 
		{
		if (tree.GetNInternals() == 1)//starTreeBefore == nInternalsBefore == 1
			return edgeMove.get();
		return lsm.get(); 
		}
	return treeScalerMove.get();

	//  double u = r->Uniform();
		//  unsigned sz = vProbTopoMove.size();
		//  double cum = 0.0;
		//  for (unsigned i = 0; i < sz; ++i)
		//  	{
		//  	cum += vProbTopoMove[i];
		//  	if (u <= cum)
		//  		{
		//  		return i;
		//  		}
		//  	}
		//  return sz;
	}
