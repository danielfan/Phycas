#ifndef	PHO_TREE_NODE_ITERATOR_TYPES_H
#define	PHO_TREE_NODE_ITERATOR_TYPES_H

#include "pyphy/src/ncl/misc/nxs_copy.hpp" // current location of Int2Type

class TreeNode;
	
enum PhoNdIteratorType
	{
		kPreorder = 0,
		kPostorder,
		kTips,
		kAwayFromNode,
		kTowardNode
	};
	
template <class T, typename U> class pho_nd_iterator;

typedef pho_nd_iterator<TreeNode, 		Int2Type<kPreorder> > 		pre_iterator;
typedef pho_nd_iterator<const TreeNode, Int2Type<kPreorder> > 		const_pre_iterator;
typedef pho_nd_iterator<TreeNode, 		Int2Type<kPostorder> > 		post_iterator;
typedef pho_nd_iterator<const TreeNode, Int2Type<kPostorder> > 		const_post_iterator;
typedef pho_nd_iterator<TreeNode, 		Int2Type<kTips> > 			tip_iterator;
typedef pho_nd_iterator<const TreeNode, Int2Type<kTips> > 			const_tip_iterator;
typedef pho_nd_iterator<TreeNode, 		Int2Type<kAwayFromNode> > 	away_iterator;
typedef pho_nd_iterator<const TreeNode, Int2Type<kAwayFromNode> > 	const_away_iterator;
typedef pho_nd_iterator<TreeNode, 		Int2Type<kTowardNode> > 	to_iterator;
typedef pho_nd_iterator<const TreeNode, Int2Type<kTowardNode> > 	const_to_iterator;

#endif

