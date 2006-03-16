#ifndef PYPHY_TREE_ITERATOR_HPP
#define PYPHY_TREE_ITERATOR_HPP

#include <boost/iterator/iterator_facade.hpp>
#include "pyphy/phylogeny/basic_tree_node.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Uses the boost library's iterator_facade to construct a valid forward iterator that visits TreeNode objects in the
|	preorder sequence. This allows constructions such as
|>
|	for (preorder_iterator nd = tree.begin(); nd != tree.end(); ++nd)
|		TakesTreeNodeRef(*nd);
|>
|	See the boost documentation on iterator_facade for a nice tutorial.
*/
class preorder_iterator 
  : public boost::iterator_facade<
		preorder_iterator,				// Derived
		TreeNode,						// Value
		boost::forward_traversal_tag>	// CategoryOrTraversal
	{
	public:

						preorder_iterator();
		explicit		preorder_iterator(TreeNode * p);

	private:
		
		friend class boost::iterator_core_access;

		void			increment();
		bool			equal(preorder_iterator const & other) const;
	    TreeNode &		dereference() const;

	private:

		TreeNode *		nd;
	};


/*----------------------------------------------------------------------------------------------------------------------
|	Default constructor needed to satisfy the forward traversal iterator requirements. Initializes data member `nd' to 
|	0.
*/
inline preorder_iterator::preorder_iterator() : nd(0)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor that sets data member `nd' to the supplied pointer `p'.
*/
inline preorder_iterator::preorder_iterator(TreeNode * p) : nd(p)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides one of the core operations required for forward iterators.
*/
inline void preorder_iterator::increment()
	{
	nd = nd->GetNextPreorder();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides one of the core operations required for forward iterators.
*/
inline bool preorder_iterator::equal(preorder_iterator const & other) const
    {
    return this->nd == other.nd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Provides one of the core operations required for forward iterators.
*/
inline TreeNode & preorder_iterator::dereference() const
    {
    return *nd;
    }

}	// namespace phycas

#endif

