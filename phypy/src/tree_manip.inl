/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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
	PHYCAS_ASSERT(t);
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
