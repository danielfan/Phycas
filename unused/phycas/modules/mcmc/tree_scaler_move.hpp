#ifndef PHO_TREESCALERMOVE_H
#define PHO_TREESCALERMOVE_H

#include "phycas/trees/tree_move.hpp"

class TreeNode;
class Tree;

/*----------------------------------------------------------------------------------------------------------------------
|	A TreeScalerMove changes the length of all edges in the tree. Each edge is set to the value x = s*x0, where x0 is
|	the edge's original length and s is a single scaling factor that applies to all edges and is calculated as follows:
|>
|	s = exp(`lambda'*(u - 0.5))
|<
|	where u is a Uniform(0,1) random deviate. The Hastings ratio for this move is s^{n}, where n is the number of edge
|	lengths in the tree.
*/
class TreeScalerMove : public TreeMove 
	{
	public:
							TreeScalerMove(Tree *t);
		virtual				~TreeScalerMove();

		void				SetLambda(double x);
		double				GetLambda() const;

		// These are pure virtual functions in the MCMCMove base class
		//
		double				GetLnHastingsRatio() const;
		double				GetLnJacobian() const;
		void				ProposeNewState();
		void				Revert();
		void				Accept();

	private:
		double				lambda;			/* larger values result in changes of greater magnitude: CV = sqrt[(lambda/2) - 1] */
		double				scalingFactor;	/* scaling factor saved in case revert is necessary */

	private:
		TreeScalerMove		&operator=(const TreeScalerMove &);	// never use - don't define
	};
typedef boost::shared_ptr<TreeScalerMove> TreeScalerMoveShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member `lambda', which is the tuning parameter for this move.
*/
inline void TreeScalerMove::SetLambda(
  double x)	/**< is the new value for the tuning parameter `lambda' */
	{
	lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the data member `lambda', which is the tuning parameter for this move.
*/
inline double TreeScalerMove::GetLambda() const
	{
	return lambda;
	}

#endif
