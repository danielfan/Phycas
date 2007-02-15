/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef	PHO_TREE_NODE_ITERATOR_TYPES_H
#define	PHO_TREE_NODE_ITERATOR_TYPES_H

#include "phycas/src/ncl/misc/nxs_copy.hpp" // current location of Int2Type

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

