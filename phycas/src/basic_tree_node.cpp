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

#include <cfloat>
#include "phycas/src/basic_tree_node.hpp"
#include "phycas/src/phycas_string.hpp"

using namespace phycas;

const double TreeNode::edgeLenEpsilon = 1.e-8;
const double TreeNode::edgeLenDefault = 0.1;
const double TreeNode::edgeLenInitValue = DBL_MAX;
const unsigned TreeNode::nodeNumInitValue = UINT_MAX;

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of immediate descendants (should be 0 for a tip node, 1 for the root (tip) node, and 2 or more for
|	internal nodes.
*/
unsigned TreeNode::CountChildren() const
	{
	unsigned nDescendants = 0;
	for (const TreeNode * child = GetLeftChildConst(); child != NULL; child = child->GetRightSibConst())
		 ++nDescendants;
	return nDescendants;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Not yet documented.
*/
void TreeNode::AppendNodeInfo(std::string &s, bool num_and_name_only)
	{
	std::string tmpstr;
	const std::string &nm = GetNodeName();

	if (nodeNum == TreeNode::nodeNumInitValue)
		{
		tmpstr << "<no number>";
		}
	else
		{
		tmpstr = str(boost::format("%d") % nodeNum);
		}

	if (num_and_name_only)
		{
		if (IsTipRoot())
			{
			s << "{";
			s << tmpstr;
			s << "}";
			}
		else if (IsTip())
			{
			s << "(";
			s << tmpstr;
			s << ")";
			}
		else if (IsInternal())
			{
			s << "[";
			s << tmpstr;
			s << "]";
			}
		else
			{
			s << "<***error***";
			s << tmpstr;
			s << "***error***>";
			}

		if (nm.length() > 0)
			{
			s << " \'";
			s << nm;
			s << "\'";
			}
		}
	else
		{
		// Node number
		s << "\nNode number: ";
		s << tmpstr;

		// Node status
		s << "\nNode status: ";
		if (IsTip())
			s << "tip";
		else if (IsTipRoot())
			s << "root";
		else if (IsInternal())
			s << "internal";
		else
			s << "unknown (this is a bug)";

		// Node name
		s << "\nNode name: ";
		if (nm.length() > 0)
			s << nm;
		else
			s << "<unnamed>";

		// lChild
		s << "\nLeft child: ";
		if (GetLeftChild() != NULL)
			lChild->AppendNodeInfo(s, true);
		else
			s << "<no left child>";

		// rSib
		s << "\nRight sibling: ";
		if (GetRightSib() != NULL)
			rSib->AppendNodeInfo(s, true);
		else
			s << "<no right sib>";

		// par
		s << "\nParent: ";
		if (GetParent() != NULL)
			par->AppendNodeInfo(s, true);
		else
			s << "<no parent>";

		// nextPreorder
		s << "\nNext preordert: ";
		if (GetNextPreorder() != NULL)
			nextPreorder->AppendNodeInfo(s, true);
		else
			s << "<no next preorder>";

		// prevPreorder
		s << "\nNext postorder: ";
		if (GetNextPostorder() != NULL)
			prevPreorder->AppendNodeInfo(s, true);
		else
			s << "<no next postorder>";

		// Edge length
		s << "\nEdge length: ";
		if (edgeLen != edgeLenInitValue)
			{
			tmpstr = str(boost::format("%f") % GetEdgeLen());
			s << tmpstr;
			}
		else
			s << "<no edge length>";
		}
	}
