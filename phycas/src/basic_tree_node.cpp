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
#include "phycas/src/basic_tree.hpp"

using namespace phycas;

#if POLPY_NEWWAY
const double TreeNode::edgeLenEpsilon = 1.e-12;  
#else
const double TreeNode::edgeLenEpsilon = 1.e-8;
#endif
const double TreeNode::edgeLenDefault = 0.1;
const double TreeNode::edgeLenInitValue = DBL_MAX;
const unsigned TreeNode::nodeNumInitValue = UINT_MAX;

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the edge length (value of `edgeLen' data member multiplied by the tree's scaler value).
*/
double TreeNode::GetEdgeLen() const
	{
    double scale = tree->GetTreeScale();
	return edgeLen*scale;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the split associated with this node.
*/
Split & TreeNode::GetSplit()
	{
	return split;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows write access to protected data member `edgeLen'.
*/
void TreeNode::SetEdgeLen(
  double x)							/**< is the new edge length value */
	{
    // The outside world doesn't know about the tree scaling factor, so the value x
    // will be an unscaled edge length
    double scale = tree->GetTreeScale();
    PHYCAS_ASSERT(scale > 0.0);
    double newlen = x/scale;
	edgeLen = (newlen < TreeNode::edgeLenEpsilon ? TreeNode::edgeLenEpsilon : newlen);

#if 0	// if reinstated, also reinstate code in TreeLikelihood::calcLnL
	if (IsInternal())
		{
		SelectNode();
		}
	else
		{
		// Internal nodes are the only ones that matter during the likelihood calculation,
		// so make sure at least one internal node gets selected. It is possible that this
		// tip node is serving as the root node, so need to make sure it has a parent before
		// trying to select the parent.
		if (par)
			par->SelectNode();
		}
#endif
	}

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
void TreeNode::AppendNodeInfo(std::string &s, bool num_and_name_only) const
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
		if (GetLeftChildConst() != NULL)
			lChild->AppendNodeInfo(s, true);
		else
			s << "<no left child>";

		// rSib
		s << "\nRight sibling: ";
		if (GetRightSibConst() != NULL)
			rSib->AppendNodeInfo(s, true);
		else
			s << "<no right sib>";

		// par
		s << "\nParent: ";
		if (GetParentConst() != NULL)
			par->AppendNodeInfo(s, true);
		else
			s << "<no parent>";

		// nextPreorder
		s << "\nNext preordert: ";
		if (GetNextPreorderConst() != NULL)
			nextPreorder->AppendNodeInfo(s, true);
		else
			s << "<no next preorder>";

		// prevPreorder
		s << "\nNext postorder: ";
		if (GetNextPostorderConst() != NULL)
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

std::string	TreeNode::briefDebugReport(
  unsigned verbosity) const
	{
	std::string nm = GetNodeName();
	std::string tmpstr;
    unsigned namelen = (unsigned)nm.length();

    if (verbosity == 0)
        {
        if (namelen == 0)
            {
		    if (IsTip())
                tmpstr = str(boost::format("(%d)") % GetNodeNumber());
            else
                tmpstr = str(boost::format("[%d]") % GetNodeNumber());
            }
        }
    else    // verbosity > 0
        {
        //     a (0) -> ? [5] -> b (1) -> ? [7] -> c (2) -> ? [6] -> d (3) -> e (4)
        if (namelen == 0)
            nm = "?";
		if (IsTip())
            tmpstr = str(boost::format(" (%d)") % GetNodeNumber());
        else
            tmpstr = str(boost::format(" [%d]") % GetNodeNumber());
        }

	return str(boost::format("%s%s") % nm % tmpstr);
	}

std::string	TreeNode::oneLineDebugReport() const
	{
	const std::string self_str = this->briefDebugReport();
	std::string lch_str;
	if (lChild)
		lch_str = lChild->briefDebugReport();
	else
		lch_str = "NULL";
	std::string par_str;
	if (par)
		par_str = par->briefDebugReport();
	else
		par_str = "NULL";
	std::string sib_str;
	if (rSib)
		sib_str = rSib->briefDebugReport();
	else
		sib_str = "NULL";
	std::string n_str;
	if (nextPreorder)
		n_str = nextPreorder->briefDebugReport();
	else
		n_str = "NULL";
	std::string prev_str;
	if (prevPreorder)
		prev_str = prevPreorder->briefDebugReport();
	else
		prev_str = "NULL";
	return str(boost::format("%s | lchild=%s | par=%s | rsib=%s | next=%s | prev=%s") % self_str % lch_str % par_str % sib_str % n_str % prev_str);
	}
