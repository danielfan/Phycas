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
class EdgeEndpoints : public std::pair<TreeNode *, TreeNode *>
	{
	public:
		EdgeEndpoints(): std::pair<TreeNode *, TreeNode *>(NULL, NULL)
			{
			}
			
		EdgeEndpoints(TreeNode * focal_nd, TreeNode * focal_neighbor)
			:std::pair<TreeNode *, TreeNode *>(focal_nd, focal_neighbor)
			{
			assert(focal_nd == NULL || focal_neighbor == NULL || (focal_neighbor->GetParent() == focal_nd || focal_nd->GetParent() == focal_neighbor));
			}
		TreeNode * getFocalNode() const
			{
			return first;
			}
		TreeNode * getFocalNeighbor() const
			{
			return second;
			}
		void setFocalNeighbor(TreeNode * n)
			{
			second = n;
			assert(n == NULL || (n->GetParent() == first || first->GetParent() == n));
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
		
		double getEdgeLength() const
			{
			TreeNode * c = this->getActualChild();
			if (c == NULL)
				{
				assert(c != NULL);
				}
			return c->GetEdgeLen();
			}
		
		operator bool () const
			{
			return first != 0;
			}
	};


}	// namespace phycas

#endif

