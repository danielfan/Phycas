#ifndef PYPHY_EDGE_ENDPOINTS_HPP
#define PYPHY_EDGE_ENDPOINTS_HPP

#include "pyphy/phylogeny/basic_tree_node.hpp"

namespace phycas
{

//typedef std::pair<TreeNode*, TreeNode*> EdgeEndpoints;
/*----------------------------------------------------------------------------------------------------------------------
|	EdgeEndpoints is a pair of adjacent nodes in the tree. 
|	For the purpose of a calculation (e.g. one that works up the tree), an edge in the tree maybe inverted wrt its 
| 		internal data structure (a child node maybe the "effective parent" of its actual parent node).
*/
class EdgeEndpoints
	{
	public:
		TreeNode * first; // effectiveChild (name legacy of when EdgeEndpoints was a pair)
		TreeNode * second; // effectiveParent (name legacy of when EdgeEndpoints was a pair)

		EdgeEndpoints()
			:first(NULL),
			second(NULL)
			{
			}
			
		EdgeEndpoints(TreeNode * effectiveChild, TreeNode * effectiveParent)
			:first(effectiveChild),
			second(effectiveParent)
			{
			assert(!(effectiveParent == NULL ^ effectiveChild == NULL));
			assert(effectiveParent == NULL || (effectiveChild->GetParent() == effectiveParent || effectiveParent->GetParent() == effectiveChild));
			}
		TreeNode * getEffectiveChild() const
			{
			return first;
			}
		TreeNode * getEffectiveParent() const
			{
			return second;
			}
		TreeNode * getActualChild() const
			{
			if (first == NULL)
				return NULL;
			if (first->GetParent() == second)
				return first;
			assert(second->GetParent() == first);
			return second;
			}
			
		TreeNode * getActualParent() const
			{
			if (first == NULL)
				return NULL;
			if (first->GetParent() == second)
				return second;
			assert(second->GetParent() == first);
			return first;
			}
		
		 operator bool () const
			{
			return first != 0;
			}
	};


}	// namespace phycas

#endif

