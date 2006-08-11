#ifndef PHO_LARGETSIMONMOVE_H
#define PHO_LARGETSIMONMOVE_H

#include <cmath>
#include "phycas/trees/tree_move.hpp"

class TreeNode;
class Tree;

class LargetSimonMove : public TreeMove
	{
	private:
		double		lambda;	/* factor used in modifying backbone length */

		TreeNode	*ndX;	/* node at one end of segment involved in move; used by Revert to undo a move */
		TreeNode	*ndY;	/* one of two node in middle of segment involved in move; used by Revert to undo a move */
		TreeNode	*ndZ;	/* node at other end (from ndX) of segment involved in move; used by Revert to undo a move */
		double		origX;	/* original length of ndX's branch; used by Revert to undo a move */
		double		origY;	/* original length of ndX's branch; used by Revert to undo a move */
		double		origZ;	/* original length of ndX's branch; used by Revert to undo a move */
		TreeNode	*swap1;	/* first of the two nodes involved in an NNI swap; NULL if no swap was performed; used by Revert to undo a move */
		TreeNode	*swap2;	/* second of the two nodes involved in an NNI swap; NULL if no swap was performed; used by Revert to undo a move */

		double		m;		/* original 3-segment length; needed for computing Hastings ratio */
		double		mstar;	/* modified 3-segment length; needed for computing Hastings ratio */

		bool		topologyChanged;	/* if true, last proposal changed topology */
		unsigned	whichcase;			/* which of the eight possible cases was tried last */

		void		Reset();

	public:
						LargetSimonMove(Tree *,double lam = 1.0);
						virtual ~LargetSimonMove() {}

		unsigned		GetWhichCase() const;
		void			SetLambda(double x);
		double			GetLambda() const;
		bool			TopologyChanged() const;

		// These are pure virtual functions in the MCMCMove base class
		//
		double		GetLnHastingsRatio() const;
		double		GetLnJacobian() const;
		void		ProposeNewState();
		void		Revert();
		void		Accept();

	private :
		LargetSimonMove &operator=(const LargetSimonMove &);//never use	- don't define
	};
typedef boost::shared_ptr<LargetSimonMove> LargetSimonMoveShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Inlined function that returns the current value of private variable whichcase, which is the index of whichever of
|	the eight possible cases was used last.
*/
inline unsigned LargetSimonMove::GetWhichCase() const
	{
	return whichcase;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inlined function that returns the natural log of the Hastings ratio for this move. The Hastings ratio is (mstar/m)^3
|	where mstar is the new length of the modified three-edge segment and m is the length before the move is proposed. 
|	Assumes that the tree is not a star tree.
*/
inline double LargetSimonMove::GetLnHastingsRatio() const
	{
#if defined(NCM_HACK)
	return 0.0;
#else
	NXS_ASSERT(tree->GetNInternals() > 1);
	return 3.0 * (std::log(mstar/m));
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inlined function that returns the natural log of the Jacobian for this move.
*/
inline double LargetSimonMove::GetLnJacobian() const
	{
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inlined function that provides read access to the private data member lambda.
*/
inline double LargetSimonMove::GetLambda() const
	{
	return lambda;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inlined function that provides read access to the private data member topologyChanged.
*/
inline bool	LargetSimonMove::TopologyChanged() const
	{
	return topologyChanged;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inlined function that provides a way of setting the value for the private data member lambda.
*/
inline void LargetSimonMove::SetLambda(double x)
	{
	lambda = x;
	}

#endif
