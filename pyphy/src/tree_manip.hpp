#ifndef TREE_MANIP_HPP
#define TREE_MANIP_HPP

//POL This file is a modified version of phycasdev/phycas/tree/tree_manip.hpp
// I am gradually incorporating the functions needed by the new Phycas

#include <cassert>
#include "boost/shared_ptr.hpp"

class ProbabilityDistribution;
typedef boost::shared_ptr<ProbabilityDistribution> ProbDistShPtr;

namespace phycas
{

class TreeNode;

class Tree;
typedef boost::shared_ptr<Tree> TreeShPtr;

class TreeID;
typedef boost::shared_ptr<TreeID> TreeIDShPtr;

class Lot;
typedef boost::shared_ptr<Lot> LotShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	The TreeManip class manipulates Tree objects (e.g. topological rearrangements needed for tree searching or MCMC).
*/
class TreeManip
	{
	public:
					TreeManip();
					TreeManip(TreeShPtr t);

		void		setTree(TreeShPtr t);

		void		starTree(unsigned ntips, ProbDistShPtr edge_len_dist);
		void		randomTree(unsigned ntips, LotShPtr rng, ProbDistShPtr edge_len_dist, bool yule);

		void		setRandomEdgeLens(ProbDistShPtr d);

		//void		BuildTreeFromID(NxsTaxaManager *tm,const TreeID& tree_id, unsigned root_at = 0);
		//void		SimpleBuildTreeFromID(unsigned nlvs, const TreeIDShPtr tree_id, unsigned root_at = 0);

		void		NNISwap(TreeNode * swap1, TreeNode * swap2);
		void		NNISwapSpecial(TreeNode * swap1);

		enum InsertMode		/**< Determines whether nodes are inserted on right or left */
			{
			kOnRight,	/**< Insert new node as rightmost child (or as right sibling of targetSib, if targetSib specified) */
			kOnLeft		/**< Insert new node as lstmost child (or as left sibling of targetSib, if targetSib specified) */
			};

		TreeNode *	FindLeftSib(TreeNode * start);
		TreeNode *	FindRightmostChild(TreeNode * start);
		TreeNode *	FindLastPreorderInClade(TreeNode * start);

		void		DetachSubtree(TreeNode * s);
		void		InsertSubtree(TreeNode * s, TreeNode * u, InsertMode m, TreeNode * targetSib = NULL);

		void		DeleteLeaf(TreeNode * u);

		void		LChildToLSib(TreeNode * u, TreeNode * w);
		void		RChildToRSib(TreeNode * u, TreeNode * w);

		void		SibToChild(TreeNode * u, TreeNode * s, InsertMode m, TreeNode * targetSib = NULL);

	protected:

		TreeShPtr	tree;	/**< The tree to be manipulated */

	};

} // namespace phycas

#include "pyphy/src/tree_manip.inl"
	
#endif
