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

#if !defined (PHO_TREE_ITERATORS_H)
#include "phycas/src/oldphycas/node_iterator.hpp"
#include "phycas/src/oldphycas/tree.hpp"

				
inline pre_iterator Tree::pre_begin() 
	{
	return pre_iterator(firstPreorder);
	}
	
inline const_pre_iterator Tree::pre_begin() const 
	{
	return const_pre_iterator(firstPreorder);
	}
	
inline post_iterator Tree::post_begin() 
	{
	return post_iterator(lastPreorder);
	}
	
inline const_post_iterator Tree::post_begin() const 
	{
	return const_post_iterator(lastPreorder);
	}
	
inline pre_iterator Tree::pre_end() 
	{
	return pre_iterator(NULL);
	}
	
inline const_pre_iterator Tree::pre_end() const 
	{
	return const_pre_iterator(NULL);
	}
	
inline post_iterator Tree::post_end() 
	{
	return post_iterator(NULL);
	}
	
inline const_post_iterator	Tree::post_end() const 
	{
	return const_post_iterator(NULL);
	}

inline tip_iterator Tree::tips_begin() 
	{
	return tip_iterator(firstPreorder);
	}
	
inline tip_iterator Tree::tips_end() 
	{
	return tip_iterator(NULL);
	}

inline tip_iterator Tree::leaves_begin() 
	{
	tip_iterator r(firstPreorder);
	++r;
	return r;
	}
	
inline tip_iterator Tree::leaves_end() 
	{
	return tip_iterator(NULL);
	}
	
inline const_tip_iterator	Tree::tips_end() const 
	{
	return const_tip_iterator(NULL);
	}
inline const_tip_iterator	Tree::tips_begin() const 
	{
	return const_tip_iterator(firstPreorder);
	}

inline const_tip_iterator Tree::leaves_begin() const
	{
	tip_iterator r(firstPreorder);
	++r;
	return r;
	}
	
inline const_tip_iterator Tree::leaves_end() const
	{
	return tip_iterator(NULL);
	}


#endif

