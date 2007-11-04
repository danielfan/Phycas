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

#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#if defined(USING_NUMARRAY)
#	define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#endif

#include <boost/python.hpp>

//#include "phycas/force_include.h"
#if defined(POL_PHYCAS)
#	include "phycas/src/basic_lot.hpp"
#else
#	include "phycas/rand/lot.hpp"
#endif
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/tree_manip.hpp"
//#include "ncl/misc/string_extensions.hpp"
//#include "ncl/nxs_exception.hpp"
#include "phycas/src/cond_likelihood.hpp"
#include "phycas/src/cond_likelihood_storage.hpp"
#include "phycas/src/tip_data.hpp"
#include "phycas/src/internal_data.hpp"

#include "phycas/src/xphylogeny.hpp"

using namespace boost::python;
using namespace phycas;

void translateXPhylogeny(const XPhylogeny &e)
	{
    // Use the Python 'C' API to set up an exception object
    PyErr_SetString(PyExc_Exception, e.what());
    }

BOOST_PYTHON_MODULE(_Phylogeny)
{
	// These function pointers are used to distinguish const from non-const member functions with the same name
	// see http://www.boost.org/libs/python/doc/tutorial/doc/html/python/functions.html#python.overloading
	const TipData *			(TreeNode::*GetTipDataConst)() const		= &TreeNode::GetTipData;
	const InternalData *	(TreeNode::*GetInternalDataConst)()	const	= &TreeNode::GetInternalData;

	class_<TreeNode>("TreeNodeBase", no_init)
		.def("getX", &TreeNode::GetX)
		.def("getY", &TreeNode::GetY)
		.def("setX", &TreeNode::SetX)
		.def("setY", &TreeNode::SetY)
		.def("isSelected", &TreeNode::IsSelected)
		.def("selectNode", &TreeNode::SelectNode)
		.def("unselectNode", &TreeNode::UnselectNode)
		.def("isRoot", &TreeNode::IsTipRoot)
		.def("isTip", &TreeNode::IsTip)
		.def("isInternal", &TreeNode::IsInternal)
		.def("getLeftChild", &TreeNode::GetLeftChild, return_internal_reference<>())
		.def("getRightSib", &TreeNode::GetRightSib, return_internal_reference<>())
		.def("getParent", &TreeNode::GetParent, return_internal_reference<>())
		.def("getNodeNumber", &TreeNode::GetNodeNumber)
		.def("getNodeName", &TreeNode::GetNodeName, return_value_policy<copy_const_reference>())
		.def("getEdgeLen", &TreeNode::GetEdgeLen)
		.def("setEdgeLen", &TreeNode::SetEdgeLen)
		.def("getNextPreorder", &TreeNode::GetNextPreorder, return_internal_reference<>())
		.def("getNextPostorder", &TreeNode::GetNextPostorder, return_internal_reference<>())
		.def("getTipData", GetTipDataConst, return_internal_reference<>())
		.def("getInternalData", GetInternalDataConst, return_internal_reference<>())
		;

	class_<Tree, boost::shared_ptr<Tree> >("TreeBase")
		.def("getNNodes", &Tree::GetNNodes)
		.def("getNTips", &Tree::GetNTips)
		.def("getNObservables", &Tree::GetNObservables)
		.def("getNInternals", &Tree::GetNInternals)
		.def("getFirstPreorder", &Tree::GetFirstPreorder, return_internal_reference<>())
		.def("getFirstPostorder", &Tree::GetLastPreorder, return_internal_reference<>())
		.def("isRooted", &Tree::IsRooted)
		.def("hasEdgeLens", &Tree::HasEdgeLens)
		.def("clear", &Tree::Clear)
		.def("buildFromString", &Tree::BuildFromString)
		.def("edgeLenSum", &Tree::EdgeLenSum)
		.def("setAllEdgeLens", &Tree::SetAllEdgeLens)
		.def("debugWalkTree", &Tree::DebugWalkTree)
		.def("debugCheckTree", &Tree::DebugCheckTree)
		.def("rerootAtTip", &Tree::RerootAtTip)
		.def("makeNewick", &Tree::MakeNewick)
		.def("tipNumbersSetUsingNames", &Tree::TipNumbersSetUsingNames)
		.def("rectifyNumbers", &Tree::RectifyNumbers)
		.def("rectifyNames", &Tree::RectifyNames)
		.def("refreshPreorder", &Tree::RefreshPreorder)
		.def("selectAllNodes", &Tree::SelectAllNodes)
		.def("unselectAllNodes", &Tree::UnselectAllNodes)
		.def("calcTotalHeight", &Tree::calcTotalHeight)
		.def("here", &Tree::DebugHere)
#if POLPY_NEWWAY
		.def("debugMode", &phycas::Tree::debugMode)
#endif
		;

	class_<TreeManip>("TreeManipBase", init<boost::shared_ptr<Tree> >())
		.def("starTree", &TreeManip::starTree)
		.def("randomTree", &TreeManip::randomTree)
		.def("setRandomEdgeLens", &TreeManip::setRandomEdgeLens)
		;

	register_exception_translator<XPhylogeny>(&translateXPhylogeny);
}
