#ifndef PYPHY_EDGE_ENDPOINTS_HPP
#define PYPHY_EDGE_ENDPOINTS_HPP

#include "pyphy/src/basic_tree_node.hpp"

namespace phycas
{

//typedef std::pair<TreeNode*, TreeNode*> EdgeEndpoints;
/*----------------------------------------------------------------------------------------------------------------------
|	EdgeEndpoints is a pair of adjacent nodes in the tree. 
|	For the purpose of a calculation (e.g. one that works up the tree), an edge in the tree maybe inverted wrt its 
| 		internal data structure (a child node maybe the "effective parent" of its actual parent node).
*/
template<typename T>
class GenericEdgeEndpoints : public std::pair<T, T>
	{
	public:
		GenericEdgeEndpoints(): std::pair<T, T>(NULL, NULL)
			{
			}
			
		GenericEdgeEndpoints(T focal_nd, T focal_neighbor)
			:std::pair<T, T>(focal_nd, focal_neighbor)
			{
			PHYCAS_ASSERT(focal_nd == NULL || focal_neighbor == NULL || (focal_neighbor->GetParentConst() == focal_nd || focal_nd->GetParentConst() == focal_neighbor));
			}
		T  getFocalNode() const
			{
			return this->first;
			}
		T  getFocalNeighbor() const
			{
			return this->second;
			}
		void setFocalNeighbor(T n)
			{
			this->second = n;
			PHYCAS_ASSERT(n == NULL || (n->GetParentConst() == this->first || this->first->GetParentConst() == n));
			}

		T 	getActualChild() const
			{
			if (this->first == NULL)
				return NULL;
			if (this->first->GetParentConst() == this->second)
				return this->first;
			PHYCAS_ASSERT(this->second->GetParentConst() == this->first);
			return this->second;
			}
			
		T 	getActualParent() const
			{
			if (this->first == NULL)
				return NULL;
			if (this->first->GetParentConst() == this->second)
				return this->second;
			PHYCAS_ASSERT(this->second->GetParentConst() == this->first);
			return this->first;
			}
		
		double getEdgeLength() const
			{
			T c = this->getActualChild();
			if (c == NULL)
				{
				PHYCAS_ASSERT(c != NULL);
				}
			return c->GetEdgeLen();
			}
		
		operator bool () const
			{
			return this->first != 0;
			}
	};
typedef GenericEdgeEndpoints<TreeNode *> EdgeEndpoints;
typedef GenericEdgeEndpoints<const TreeNode *> ConstEdgeEndpoints;

}	// namespace phycas

#endif

