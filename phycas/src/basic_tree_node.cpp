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

const double TreeNode::edgeLenEpsilon = 1.e-12;  
const double TreeNode::edgeLenDefault = 0.1;
const double TreeNode::edgeLenInitValue = DBL_MAX;
const unsigned TreeNode::nodeNumInitValue = UINT_MAX;

//static unsigned debug_treenode_number = 0;


/*----------------------------------------------------------------------------------------------------------------------
|	Only sets par pointer and either lChild (of parent) or rSib (of left sib).
*/
TreeNode * TreeNode::AddChild(TreeNode *nd)
    {
    TreeNode * leftSib = FindRightmostChild();
    if (leftSib == 0L)
        lChild = nd;
    else
        leftSib->rSib = nd;
    nd->par = this;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes all pointers to 0, nodeNum to TreeNode::nodeNumInitValue, and edgeLen to TreeNode::edgeLenInitValue.
*/
TreeNode::TreeNode() 
  :nodeNum(TreeNode::nodeNumInitValue),
  edgeLen(TreeNode::edgeLenInitValue),
  lChild(0),
  par(0),
  rSib(0),
  nextPreorder(0),
  prevPreorder(0),
  //observable(false),
  support(0.0),
  tmp(0.0),
  x(0.0),
  y(0.0),
  selected(false),
  tipData(0),
  tipDataDeleter(0),
  internalData(0),
  internalDataDeleter(0)
  	{
  	//std::cerr << "=====> creating TreeNode " << (++debug_treenode_number) << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor for TreeNode. 
*/
TreeNode::~TreeNode()
	{
	//std::cerr << "In node destructor" << std::endl;
	Clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes all pointers to 0, nodeNum to TreeNode::nodeNumInitValue and edgeLen to TreeNode::edgeLenInitValue.
|	Also deletes any structures assigned to `tipData' or `internalData' using the callbacks provided when these
|	structures were allocated. Basically, returns node to its just-constructed state.
*/
void TreeNode::Clear()
	{
	//std::cerr << "tree use count = " << tree.use_count() << std::endl;
	//tree.reset();
	ResetInternalData();
	ResetTipData();

	lChild			= 0; 
	par				= 0;
	rSib			= 0;
	nextPreorder	= 0;
	prevPreorder	= 0;
	nodeNum			= TreeNode::nodeNumInitValue;
	edgeLen			= TreeNode::edgeLenInitValue;
	nodeName		= "";
	//observable		= false;
	support			= 0.0;
	tmp				= 0.0;
	x				= 0.0;
	y				= 0.0;
	selected		= false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the edge length (value of `edgeLen' data member multiplied by the tree's scaler value).
*/
double TreeNode::GetEdgeLen() const
	{
	return edgeLen;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the split associated with this node.
*/
const Split & TreeNode::GetSplitConst() const
	{
	return split;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the split associated with this node.
*/
Split & TreeNode::GetSplit()
	{
	return split;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `edgeLen' to the product of the current value of `edgeLen' and the supplied `scaling_factor'.
*/
void TreeNode::ScaleEdgeLen(
  double scaling_factor)			/**< is the value by which the edge length will be multiplied */
	{
    PHYCAS_ASSERT(scaling_factor > 0.0);
    double new_edgeLen = edgeLen*scaling_factor;
    SetEdgeLen(new_edgeLen);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Allows write access to protected data member `edgeLen'.
*/
void TreeNode::SetEdgeLen(
  double x)							/**< is the new edge length value */
	{
	edgeLen = (x < TreeNode::edgeLenEpsilon ? TreeNode::edgeLenEpsilon : x);
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
