/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

#ifndef PHO_TREEMANIP_H
#define PHO_TREEMANIP_H

class TreeNode;
class Tree;

#if defined(POL_PHYCAS)
namespace phycas
	{
	class Lot;
	}
#else
class Lot;
#endif

class TreeID;
class NxsTaxaManager;	// 22jun2005

class TreeManip
	{
	public:
					TreeManip(Tree *t);

#if defined(POL_PHYCAS)
		void		SimpleBuildStarTree(unsigned nlvs, unsigned root_at = 0, phycas::Lot * r = NULL, double lambda = 1.0);	//POL added root_at 22-Oct-2004
		void		SimpleBuildYuleTree(phycas::Lot & r, double lambda, unsigned nlvs);
#else
		void		SimpleBuildStarTree(unsigned nlvs, unsigned root_at = 0, Lot *r = NULL, double lambda = 1.0);	//POL added root_at 22-Oct-2004
		void		SimpleBuildYuleTree(Lot &r, double lambda, unsigned nlvs);
#endif
		void		BuildTreeFromID(NxsTaxaManager *tm,const TreeID& tree_id, unsigned root_at = 0);	//POL added 22jun2005
		void		SimpleBuildTreeFromID(unsigned nlvs,const TreeID& tree_id, unsigned root_at = 0);	//POL added root_at 22-Oct-2004

		enum InsertMode		/* Determines whether nodes are inserted on right or left */
			{
			kOnRight,	/* insert new node as rightmost child (or as right sibling of targetSib, if targetSib specified) */
			kOnLeft		/* insert new node as lstmost child (or as left sibling of targetSib, if targetSib specified) */
			};

	protected:
		Tree		*tree;

		void		DetachSubtree(TreeNode *s);
		void		InsertSubtree(TreeNode *s, TreeNode *u, InsertMode m, TreeNode *targetSib = NULL);

		void		DeleteLeaf(TreeNode *u);

		void		NNISwap(TreeNode *swap1, TreeNode *swap2);
		void		NNISwapSpecial(TreeNode *swap1);

		void		LChildToLSib(TreeNode *u, TreeNode *w);
		void		RChildToRSib(TreeNode *u, TreeNode *w);

		void		SibToChild(TreeNode *u, TreeNode *s, InsertMode m, TreeNode *targetSib = NULL);
	};
	
inline TreeManip::TreeManip(
  Tree *t)
  	{
	PHYCAS_ASSERT(t != NULL);
	tree = t;
  	}


#endif
