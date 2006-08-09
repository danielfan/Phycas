#ifndef PYPHY_TOWARD_NODE_ITERATOR_HPP
#define PYPHY_TOWARD_NODE_ITERATOR_HPP

#include <utility>
#include <stack>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/function.hpp>
#include "pyphy/phylogeny/basic_tree_node.hpp"
#include "pyphy/phylogeny/edge_endpoints.hpp"

namespace phycas
{
typedef boost::function< bool (TreeNode *, TreeNode *) > NodeValidityChecker;
		
/*----------------------------------------------------------------------------------------------------------------------
|	Iterates through EdgeEndpoints. The ordering of nodes returned is "post-order" (post-order if the node used to 
|	construct the iterator were the root), but the creation of the iterator involves a pre-order traversal to allow 
|	subtrees to be omitted if they are valid according to the NodeValidityChecker functor.
*/
class effective_postorder_edge_iterator 
  : public boost::iterator_facade<
		effective_postorder_edge_iterator,	// Derived
		EdgeEndpoints,						// Value
		boost::forward_traversal_tag>		// CategoryOrTraversal
	{
	public:

						effective_postorder_edge_iterator();
		explicit		effective_postorder_edge_iterator(TreeNode * effectiveRoot, NodeValidityChecker f /* = NodeValidityChecker() <-- this causes internal compiler error in VC */);

	private:
		
		friend class boost::iterator_core_access;

		void						increment();
		bool				  		equal(effective_postorder_edge_iterator const & other) const;
	    EdgeEndpoints &				dereference() const;

	private:
		void						BuildStackFromNodeAndSiblings(TreeNode *, const TreeNode *);

		NodeValidityChecker			isValidChecker;
		std::stack<EdgeEndpoints>	edge_stack;
		mutable EdgeEndpoints		topOfEdgeStack;
		TreeNode *					effectiveRoot;
	};

class EffectivePostorderEdgeContainer
	{
	public:
		EffectivePostorderEdgeContainer(TreeNode * effectiveRoot, NodeValidityChecker f)
			:beginIt(effectiveRoot, f)
			{}
			
		effective_postorder_edge_iterator	begin()
			{
			return beginIt;
			}

		effective_postorder_edge_iterator	end()
			{
			return effective_postorder_edge_iterator();
			}
	private:
		effective_postorder_edge_iterator beginIt;
	};
	
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline effective_postorder_edge_iterator::effective_postorder_edge_iterator() : effectiveRoot(NULL)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline effective_postorder_edge_iterator::effective_postorder_edge_iterator(
  TreeNode * effRoot, 		/**< "focal" node of the iteration (order of nodes will be postorder if this node were the root) */
  NodeValidityChecker f)	/**< functor that takes two Node pointers and returns true if the iteration in this subtree should be stopped) */
  : isValidChecker(f), effectiveRoot(effRoot)
	{
	assert(!isValidChecker.empty());
	assert(effectiveRoot != NULL);

	// Build that part of edge_stack based on nodes above effectiveRoot
	BuildStackFromNodeAndSiblings(effectiveRoot->GetLeftChild(), NULL);

	// Now finish the job by considering nodes below effectiveRoot
	TreeNode * currAvoidNode = effectiveRoot;
	for (TreeNode * currAnc = effectiveRoot->GetParent(); currAnc != NULL; currAnc = currAnc->GetParent())
		{
		if (!isValidChecker.empty() && isValidChecker(currAnc, currAvoidNode))
			break;

		// 
		edge_stack.push(EdgeEndpoints(currAnc, currAvoidNode));
		BuildStackFromNodeAndSiblings(currAnc->GetLeftChild(), currAvoidNode);
		currAvoidNode = currAnc;
		}

	// Make sure we look like a default-constructed object if there are no edges in edge_stack
	// This allows the iterator to determine when it has reached the end
	if (edge_stack.empty())
		effectiveRoot = NULL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline void effective_postorder_edge_iterator::BuildStackFromNodeAndSiblings(TreeNode * curr, const TreeNode * nodeToSkip)
	{
	// Handle case in which curr equals nodeToSkip
	if (nodeToSkip != NULL && curr == nodeToSkip)
		curr = curr->GetRightSib();
	if (curr == NULL)
		return;

	// Create a temporary stack of nodes to remember
	std::stack<TreeNode *> nd_stack;

	// Visit all nodes in the subtree extending up from curr's parent (with the exception 
	// of the subtree whose root is nodeToSkip)
	for (;;)
		{
		TreeNode * par = curr->GetParent();
		if ((!isValidChecker.empty()) && isValidChecker(curr, par))
			{
			// edge was valid:
			//   - let curr be curr's next sibling
			curr = curr->GetRightSib();
			if (nodeToSkip != NULL && curr == nodeToSkip)
				curr = curr->GetRightSib();
			}
		else
			{
			// edge was not valid:
			//   - push edge onto edge stack
			//   - if curr has a sibling, push that sibling onto nd_stack (i.e. remember it so we can deal with it later)
			//   - let curr be curr's left child
			edge_stack.push(EdgeEndpoints(curr, par));
			TreeNode * r = curr->GetRightSib();
			if (r != NULL && r == nodeToSkip)
				r = r->GetRightSib();
			if (r != NULL)
				nd_stack.push(r);
			curr = curr->GetLeftChild();
			}

		if (curr == NULL)
			{
			// we've come to the end of the road for this subtree
			// let curr be node on top of nd_stack
			if (nd_stack.empty())
				break;
			curr = nd_stack.top();
			nd_stack.pop();
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline void effective_postorder_edge_iterator::increment()
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
inline bool effective_postorder_edge_iterator::equal(effective_postorder_edge_iterator const & other) const
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
inline EdgeEndpoints & effective_postorder_edge_iterator::dereference() const
    {
	assert(!edge_stack.empty());
	topOfEdgeStack = edge_stack.top();
    return topOfEdgeStack;
    }

}	// namespace phycas

#endif

