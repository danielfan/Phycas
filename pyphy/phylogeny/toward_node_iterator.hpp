#ifndef PYPHY_TOWARD_NODE_ITERATOR_HPP
#define PYPHY_TOWARD_NODE_ITERATOR_HPP

#include <utility>
#include <stack>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/function.hpp>
#include "pyphy/phylogeny/basic_tree_node.hpp"
#include "pyphy/likelihood/tree_likelihood.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	EdgeEndpoints is a pair of adjacent nodes.
*/
typedef std::pair<TreeNode*, TreeNode*> EdgeEndpoints;
typedef boost::function< bool (const TreeNode *, const TreeNode *) > NodeValidityChecker;
		
/*----------------------------------------------------------------------------------------------------------------------
|	Iterates through EdgeEndpoints 
|	The ordering of nodes returned is "post-order" (post-order if the node used to construct the iterator were the root)
|	But the creation of the iterator involves a pre-order traversal to allow subtree to be 
|		omitted if they are valid.
*/
class toward_nd_iterator 
  : public boost::iterator_facade<
		toward_nd_iterator,				// Derived
		EdgeEndpoints,					// Value
		boost::forward_traversal_tag>	// CategoryOrTraversal
	{
	public:

						toward_nd_iterator();
		explicit		toward_nd_iterator(TreeNode * effectiveRoot, NodeValidityChecker f);

	private:
		
		friend class boost::iterator_core_access;

		void					increment();
		bool				  	equal(toward_nd_iterator const & other) const;
	    EdgeEndpoints	      & dereference() const;

	private:
		void					BuildStackFromSubtree(TreeNode *);

		NodeValidityChecker			isValidChecker;
		std::stack<EdgeEndpoints>	edge_stack;
		mutable EdgeEndpoints		topOfEdgeStack;
		TreeNode 				  * effectiveRoot;
	};

class TowardNodeIteratorContainer
	{
	public:
		TowardNodeIteratorContainer(TreeNode * effectiveRoot, NodeValidityChecker f)
			:beginIt(effectiveRoot, f)
			{}
			
		toward_nd_iterator	begin()
			{
			return beginIt;
			}

		toward_nd_iterator	end()
			{
			return toward_nd_iterator();
			}
	private:
		toward_nd_iterator beginIt;
	};
	
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline toward_nd_iterator::toward_nd_iterator() : effectiveRoot(NULL)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline toward_nd_iterator::toward_nd_iterator(
  TreeNode * effRoot, 		/// "focal" node of the iteration (order of nodes will be postorder if this node were the root)
  NodeValidityChecker f)	/// functor that takes two Node pointers and returns true if the iteration in this subtree should be stopped)
  : isValidChecker(f), effectiveRoot(effRoot)
	{
	assert(effectiveRoot != NULL);
	TreeNode * leftChild = effectiveRoot->GetLeftChild();
	assert(leftChild != NULL);
	BuildStackFromSubtree(leftChild);
	TreeNode * currAvoidNode = effectiveRoot;
	for (TreeNode *currAnc = effectiveRoot->GetParent(); currAnc != NULL; currAnc = currAnc->GetParent(), currAvoidNode = currAnc)
		{
		if ((!isValidChecker.empty()) && isValidChecker(currAnc, currAvoidNode))
			break;
		edge_stack.push(EdgeEndpoints(currAnc, currAvoidNode));
		for (TreeNode * c = currAnc->GetLeftChild(); c != NULL; c = c->GetLeftChild())
			{
			if (c != currAvoidNode && (isValidChecker.empty() || !isValidChecker(c, currAnc)))
				{
				edge_stack.push(EdgeEndpoints(c, currAnc));
				BuildStackFromSubtree(c);
				}
			}
		}
	if (edge_stack.empty())
		effectiveRoot = NULL;
	}


inline void toward_nd_iterator::BuildStackFromSubtree(TreeNode * curr)
	{
	std::stack<TreeNode*> nd_stack;
	for (;;)
		{
		TreeNode * par = curr->GetParent();
		if ((!isValidChecker.empty()) && isValidChecker(curr, par))
			{
			curr = curr->GetRightSib();
			if (curr == NULL)
				{
				if (nd_stack.empty())
					break;
				curr = nd_stack.top();
				nd_stack.pop();
				}
			}
		else
			{
			edge_stack.push(EdgeEndpoints(curr, par));
			TreeNode *r = curr->GetRightSib();
			if (r != NULL)
				nd_stack.push(r);
			curr = curr->GetLeftChild();
			assert(curr != NULL);
			}
		}
	}
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline void toward_nd_iterator::increment()
	{
	if (!edge_stack.empty())
		{
		edge_stack.pop();
		if (edge_stack.empty())
			effectiveRoot = NULL;
		}
	else
		{
		assert(effectiveRoot == NULL);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides one of the core operations required for forward iterators.
*/
inline bool toward_nd_iterator::equal(toward_nd_iterator const & other) const
    {
    if (this->effectiveRoot != other.effectiveRoot)
    	return false;
    if (edge_stack.empty())
    	return other.edge_stack.empty();
    return (!other.edge_stack.empty()) && (edge_stack.top().first == other.edge_stack.top().first);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Provides one of the core operations required for forward iterators.
*/
inline EdgeEndpoints & toward_nd_iterator::dereference() const
    {
	assert(!edge_stack.empty());
	topOfEdgeStack = edge_stack.top();
    return topOfEdgeStack;
    }

}	// namespace phycas

#endif

