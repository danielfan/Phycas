#if !defined(PHYCAS_TREES_TREE_INL_H)
#include "pyphy/src/oldphycas/tree.hpp"
#include "pyphy/src/oldphycas/tree_node.hpp"

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

