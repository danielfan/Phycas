#ifndef TREE_MANIP_INL
#define TREE_MANIP_INL

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	The default constructor does not set the data member `tree'. You should set `tree' using the SetTree() member
|	function before using any of the other member functions, most of which assume that `tree' points to a Tree object.
*/
inline TreeManip::TreeManip()
  	{
  	}

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor requires a single argument, which is a shared pointer to the tree to be manipulated.
*/
inline TreeManip::TreeManip(
  TreeShPtr t)	/**< Is the tree to be manipulated */
  	{
	assert(t);
	tree = t;
  	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `tree' data member to the supplied shared pointer to Tree.
*/
inline void TreeManip::setTree(TreeShPtr t)
  	{
	tree = t;
  	}

} // phycas namespace

#endif
