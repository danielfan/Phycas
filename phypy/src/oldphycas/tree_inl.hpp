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

#if !defined(PHYCAS_TREES_TREE_INL_H)
#include "phypy/src/oldphycas/tree.hpp"
#include "phypy/src/oldphycas/tree_node.hpp"

// This file contains Tree member functions that should be inlined, but which (at least on gcc3.4) require
// 	the include "tree/tree_node.hpp" statement (in previous versions of gcc tree node.hpp would have to be
// 	included in the same compilation unit, but not necessarily in the same file).  
// Rather than add this include, to tree.hpp, the functions were moved to this file.

/*----------------------------------------------------------------------------------------------------------------------
|       Calls a functor that takes a const TreeNode* for all leaves in the tree
*/
template <class FuncToCall>inline void Tree::ForAllLeaves(FuncToCall f) const
        {
        for (const TreeNode *nd = GetRoot(); nd != NULL; nd = nd->GetNextPreorder())
                {
                if (nd->IsShootTip()) ///@ is this a bug? should this be is RootOrLeaf
                        f(nd);
                }
        }

/*----------------------------------------------------------------------------------------------------------------------
|       Calls a functor that takes a TreeNode * for all leaves in the tree
*/
template <class FuncToCall>inline void Tree::ForAllLeaves(FuncToCall f)
        {
        for (TreeNode *nd = GetRoot(); nd != NULL; nd = nd->GetNextPreorder())
                {
                if (nd->IsShootTip())
                        f(nd);
                }
        }
/*----------------------------------------------------------------------------------------------------------------------
|       Traverses tree in postorder fashion, calling functor f for each node visited. This provides a way for classes 
|       such as SplitManager to have access to nodes without Tree having to know the details of the other class (which 
|       keeps Tree uncomplicated). Supplied functor should wrap a function returning a bool and taking a single TreeNode* 
|       parameter.
*/
template <class FuncToCall>
inline void Tree::PostorderTraverse(
  FuncToCall f) /**< function to be called for each node visited */
        {
        for (TreeNode *nd = GetFirstPostorder(); nd != NULL; nd = nd->GetNextPostorder())
                {
                if (!f(nd))
                        return;
                }
        }

#endif

