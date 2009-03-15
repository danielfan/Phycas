/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
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

#ifndef PYPHY_BASIC_TREE_NODE_INL
#define PYPHY_BASIC_TREE_NODE_INL

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the rSib or the parent's left child.
|	Returns NULL if the parent has only one child.
*/
inline TreeNode * TreeNode::FindNextSib()
	{
	if (rSib != NULL)
		return rSib;
	PHYCAS_ASSERT(par);
	if (par->lChild == this)
		return NULL;
	return par->lChild;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `tree' share pointer data member to the supplied value `t'. This allows subsequent access to the tree of
|   which this node is a part in order to, for example, gain access to the whole-tree scaling factor.
*/
inline void TreeNode::SetTreeShPtr(
  TreeShPtr t) /**> is a pointer to the tree (should be called just after a node is created) */
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `support'.
*/
inline float TreeNode::GetSupport()
	{
	return support;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `x'.
*/
inline float TreeNode::GetX()
	{
	return x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `observable' to true.
*/
inline void TreeNode::SetObservable()
	{
	observable = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `observable' to false.
*/
inline void TreeNode::SetUnobservable()
	{
	observable = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `support' to the value `x'.
*/
inline void TreeNode::SetSupport(
  float x)						/**< is the new support value */
	{
	support = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `x' to the value `xx'.
*/
inline void TreeNode::SetX(
  float xx)						/**< is the new x-coordinate */
	{
	x = xx;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `y'.
*/
inline float TreeNode::GetY()
	{
	return y;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `y' to the value `yy'.
*/
inline void TreeNode::SetY(
  float yy)						/**< is the new y-coordinate */
	{
	y = yy;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier that sets `nodeNum' data member to the supplied value. 
*/
inline void TreeNode::SetNodeNum(
  unsigned num)						/**< is the new node number */
	{
	nodeNum = num;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier that sets `nodeName' data member to the supplied value. 
*/
inline void TreeNode::SetNodeName(
  std::string name)					/**< is the new node name */
	{
	nodeName = name;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `lChild' data member is not zero.
*/
inline bool TreeNode::HasChildren() const
	{
	return (lChild != 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `lChild' data member equals zero.
*/
inline bool TreeNode::NoChildren() const
	{
	return (lChild == 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `par' data member is not zero.
*/
inline bool TreeNode::HasParent() const
	{
	return (par != 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `par' data member equals zero.
*/
inline bool TreeNode::NoParent() const
	{
	return (par == 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if either NoChildren() or IsTipRoot() returns true. A node is a tip node if it has degree
|	one, which is true if only one edge connects it to the rest of the tree. This is true for nodes that have no 
|	children (their only connection to the rest of the tree is through their parent), but is also true for the root node
|	(the only node in a tree that has only one child).
*/
inline bool TreeNode::IsTip() const
	{
	return NoChildren() || IsTipRoot();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if the node can potentially have associated observations. For unrooted trees, this will be
|	true for any IsTip() node, but for rooted trees, it is only true for NoChildren() nodes because the root is
|	not observable even though it is degree one. Note that a TreeNode object has no way of knowing whether the tree to
|	which it is connected is rooted or unrooted, so this function simply returns the value of the `observable' data 
|	member, and depends on the tree to set that data member correctly when the node is incorporated into the tree.
*/
inline bool TreeNode::IsObservable() const
	{
	return observable;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets flag that serves to indicate whether the node can potentially have associated observations. For unrooted trees,
|	this should be true for any IsTip() node, but for rooted trees, it should only be true for NoChildren() nodes
|	because the root is	not observable even though it is degree one. Note that a TreeNode object has no way of knowing 
|	whether the tree to which it is connected is rooted or unrooted, so this function should be used by functions that
|	build trees to set the `observable' data member correctly. One use of this data member is in deciding which nodes
|	should be equipped with structures that store a copy of a particular row in the data matrix when preparing the tree
|	for likelihood calculations.
*/
inline void TreeNode::SetObservable(bool is_observable)
	{
	observable = is_observable;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if IsTip() returns false.
*/
inline bool TreeNode::IsInternal() const
	{
	return !IsTip();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `nodeNum' equals the value `TreeNode::nodeNumInitValue', which is the value assigned to
|	the `nodeNum' data member of all newly-created nodes.
*/
inline bool TreeNode::NumberNotYetAssigned() const
	{
	return nodeNum == TreeNode::nodeNumInitValue;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `edgeLen' equals the value `TreeNode::edgeLenInitValue', which is the value assigned to
|	the `edgeLen' data member of all newly-created nodes.
*/
inline bool TreeNode::EdgeLenNotYetAssigned() const
	{
	return edgeLen == TreeNode::edgeLenInitValue;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if this is a node with no parent and no siblings. This function does not pay attention to
|	the number of children. If you know that the root node should be a tip, use TreeNode::IsTipRoot instead (that 
|	function also requires root to have just one child). If you know that the root node should be an internal node, use
|	TreeNode::IsInternalRoot instead (that function also requires that the root have at least two children).
*/
inline bool TreeNode::IsAnyRoot() const
	{
	bool parentless = !par;
	bool only_child = !rSib;
	return (parentless && only_child);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if this is a node with no parent, no siblings and only one child. Use this function if
|	you specifically want to test for the case where a tip node is serving as the root of the tree. If an internal node
|	is serving as the root of the tree, it is expected to have more than one child, and the method IsInternalRoot()
|	should be used instead.
*/
inline bool TreeNode::IsTipRoot() const	//POL: formerly IsRoot()
	{
	bool parentless = !par;
	bool only_child = !rSib;
	bool one_child = (lChild && !(lChild->rSib));
	return (parentless && only_child && one_child);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if this is a node with no parent, no siblings and at least 2 children.
*/
inline bool TreeNode::IsInternalRoot() const
	{
	bool parentless = !par;
	bool only_child = !rSib;
	bool at_least_two_children = (lChild && lChild->rSib);
	return (parentless && only_child && at_least_two_children);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the edge length (value of `nodeNum' data member).
*/
inline unsigned TreeNode::GetNodeNumber() const
	{
	return nodeNum;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the name of this node (value of `nodeName' data member).
*/
inline const std::string & TreeNode::GetNodeName() const
	{
	return nodeName;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `selected' data member. This function (as well as its cousings SelectNode and UnselectNode) may be
|	temporary. They are intended to be used for situations in which some nodes need to be selected, either in a 
|	graphical context when manipulating trees via a GUI, or in a non-graphical context, such as identification of nodes
|	for lazy recalculation. Both of these applications are probably better served by having more specific flags, 
|	however, so these functions will probably eventually be moved into other structures.
*/
inline bool TreeNode::IsSelected() const
	{
	return selected;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `selected' data member to true. (See notes under IsSelected.)
*/
inline void TreeNode::SelectNode()
	{
	selected = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `selected' data member to false. (See notes under IsSelected.)
*/
inline void TreeNode::UnselectNode()
	{
	selected = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `tipData' data member.
*/
inline TipData * TreeNode::GetTipData()
	{
	return tipData;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `internalData' data member.
*/
inline InternalData * TreeNode::GetInternalData()
	{
	return internalData;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `tipData' data member.
*/
inline const TipData * TreeNode::GetTipData() const
	{
	return tipData;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `internalData' data member.
*/
inline const InternalData * TreeNode::GetInternalData() const
	{
	return internalData;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `tipData' to its just-constructed state, namely 0. If `tipData' is non-zero, calls the associated
|	`tipDataDeleter' function to delete the object. Assumes that if `tipData' is non-zero, then `tipDataDeleter' is 
|	valid.
*/
inline void TreeNode::ResetTipData()
	{
	if (tipData)
		{
		PHYCAS_ASSERT(tipDataDeleter);
		tipDataDeleter(tipData);
		}
	tipData = 0;
	tipDataDeleter = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `internalData' to its just-constructed state, namely 0. If `internalData' is non-zero, calls the
|	associated `internalDataDeleter' function to delete the object. Assumes that if `internalData' is non-zero, then 
|	`internalDataDeleter' is valid.
*/
inline void TreeNode::ResetInternalData()
	{
	if (internalData)
		{
		PHYCAS_ASSERT(internalDataDeleter);
		internalDataDeleter(internalData);
		}
	internalData = 0;
	internalDataDeleter = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `tipData' to `d' and `tipDataDeleter' to `f'
*/
inline void TreeNode::SetTipData(
  TipData * d,						/**< is a pointer to the TipData data structure to be assigned to `tipData' */
  TreeNode::TipDataDeleter f)		/**< is the function object that knows how to delete a TipData object */
	{
	ResetTipData();
	tipData = d;
	tipDataDeleter = f;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `internalData' to `c' and `internalDataDeleter' to `f'
*/
inline void TreeNode::SetInternalData(
  InternalData * c,						/**< is a pointer to the InternalData data structure to be assigned to `internalData' */
  TreeNode::InternalDataDeleter f)		/**< is the function object that knows how to delete a InternalData object */
	{
	ResetInternalData();
	internalData = c;
	internalDataDeleter = f;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `lChild'.
*/
inline TreeNode	* TreeNode::GetLeftChild()
	{
	return lChild;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `lChild'. Useful for iterating through children of a node inside a const
|	function.
*/
inline const TreeNode * TreeNode::GetLeftChildConst() const
	{
	return lChild;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `rSib'.
*/
inline TreeNode	* TreeNode::GetRightSib()
	{
	return rSib;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `rSib'. Useful for iterating through children of a node inside a const
|	function.
*/
inline const TreeNode * TreeNode::GetRightSibConst() const
	{
	return rSib;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `par'.
*/
inline TreeNode	* TreeNode::GetParent()
	{
	return par;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `par'.
*/
inline const TreeNode	* TreeNode::GetParentConst() const
	{
	return par;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `nextPreorder'.
*/
inline TreeNode	* TreeNode::GetNextPreorder()
	{
	return nextPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `prevPreorder'.
*/
inline TreeNode	* TreeNode::GetNextPostorder()
	{
	return prevPreorder;
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `nextPreorder'.
*/
inline const TreeNode	* TreeNode::GetNextPreorderConst() const 
	{
	return nextPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `prevPreorder'.
*/
inline const TreeNode	* TreeNode::GetNextPostorderConst() const 
	{
	return prevPreorder;
	}

} // namespace phycas

#endif

