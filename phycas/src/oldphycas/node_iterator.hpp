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

#ifndef	PHO_TREE_NODE_ITERATOR_H
#define	PHO_TREE_NODE_ITERATOR_H

#include "phycas/src/oldphycas/node_iterator_types.hpp"
#include "phycas/src/oldphycas/tree_node.hpp"

template <class T, typename U> class pho_nd_iterator
	{
	public:
		explicit 	pho_nd_iterator(T *n = NULL);
		explicit 	pho_nd_iterator(T &n);
					pho_nd_iterator(const pho_nd_iterator<TreeNode, U> &n): nd(n.operator->()) 
			{}
		T *operator->() const	{	return nd;	}
		T &operator*() const	{	return *nd;	}
		pho_nd_iterator &operator++();
		pho_nd_iterator operator++(int)
			{
			pho_nd_iterator ret(*this);
			++(*this);
			return ret;
			}
		template <class ST, typename SB> bool operator==(const pho_nd_iterator<ST, SB> &r) const	{	return (r.operator->() == nd);	}
		bool operator==(const T &n) const	{	return (&n == nd); }
		bool operator==(const T *n) const	{ 	return (n == nd); }
		template <class ST, typename SB> bool operator!=(const pho_nd_iterator<ST, SB> &r) const	{	return (r.operator->() != nd);	}
		bool operator!=(const T &n) const	{	return (&n != nd); }
		bool operator!=(const T *n) const	{ 	return (n != nd); }
	private:
		T *nd;
	};

template<> inline pre_iterator &pre_iterator::operator++()
	{
	NXS_ASSERT(nd != NULL);
	nd = nd->GetNextPreorder();
	return *this;
	}
		
template<> inline const_pre_iterator &const_pre_iterator::operator++()
	{
	NXS_ASSERT(nd != NULL);
	nd = nd->GetNextPreorder();
	return *this;
	}
		
template<> inline post_iterator &post_iterator::operator++()
	{
	NXS_ASSERT(nd != NULL);
	nd = nd->GetNextPostorder();
	return *this;
	}
		
template<> inline const_post_iterator &const_post_iterator::operator++()
	{
	NXS_ASSERT(nd != NULL);
	nd = nd->GetNextPostorder();
	return *this;
	}
		
template<> inline tip_iterator &tip_iterator::operator++()
	{
	NXS_ASSERT(nd != NULL);
	nd = nd->GetNextPreorder();
	while (nd != NULL && nd->IsInternal())
		nd = nd->GetNextPreorder();
	return *this;
	}

template<> inline const_tip_iterator &const_tip_iterator::operator++()
	{
	NXS_ASSERT(nd != NULL);
	nd = nd->GetNextPreorder();
	while (nd != NULL && nd->IsInternal())
		nd = nd->GetNextPreorder();
	return *this;
	}

template <> inline tip_iterator::pho_nd_iterator(TreeNode *n)
	: nd(n) 
	{
	if (nd != NULL && nd->IsInternal())
		++(*this);
	}
	
template <> inline tip_iterator::pho_nd_iterator(TreeNode &n)
	: nd(&n) 
	{
	if (nd != NULL && nd->IsInternal())
		++(*this);
	}
	

template <> inline const_tip_iterator::pho_nd_iterator(const TreeNode *n)
	: nd(n) 
	{
	if (nd != NULL && nd->IsInternal())
		++(*this);
	}
	
template <> inline const_tip_iterator::pho_nd_iterator(const TreeNode &n)
	: nd(&n) 
	{
	if (nd != NULL && nd->IsInternal())
		++(*this);
	}
	

template <class T, typename U> inline pho_nd_iterator<T, U>::pho_nd_iterator(T *n)
	: nd(n) 
	{}
	
template <class T, typename U> inline pho_nd_iterator<T, U>::pho_nd_iterator(T &n)
    : nd(&n) 
    {}
	
		
		
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the first preorder node for this clade -- this.
*/
inline pre_iterator TreeNode::clade_pre_begin() 
	{
	return pre_iterator(this);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the first preorder node for this clade -- this.
*/
inline const_pre_iterator TreeNode::clade_pre_begin() const 
	{
	return const_pre_iterator(this);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the first postorder node for this clade -- FindLastPreorderInClade().
*/
inline post_iterator TreeNode::clade_post_begin() 
	{
	return post_iterator(FindLastPreorderInClade());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the first postorder node for this clade -- FindLastPreorderInClade().
*/
inline const_post_iterator TreeNode::clade_post_begin() const 
	{
	return const_post_iterator(FindLastPreorderInClade());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the end preorder node for this clade -- FindLastPreorderInClade()->GetNextPreorder()
*/
inline pre_iterator TreeNode::clade_pre_end() 	
	{	
	return pre_iterator(FindLastPreorderInClade()->GetNextPreorder());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the end preorder node for this clade -- FindLastPreorderInClade()->GetNextPreorder()
*/
inline const_pre_iterator TreeNode::clade_pre_end() const 
	{
	return const_pre_iterator(FindLastPreorderInClade()->GetNextPreorder());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the end postorder node for this clade -- GetNextPostorder()
*/
inline post_iterator TreeNode::clade_post_end() 
	{	
	return post_iterator(GetNextPostorder());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the end postorder node for this clade -- GetNextPostorder()
*/
inline const_post_iterator TreeNode::clade_post_end() const 
	{
	return const_post_iterator(GetNextPostorder());
	}

#endif
