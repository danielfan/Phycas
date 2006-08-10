#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle		// only needed on lewis.eeb.uconn.edu

//#include "phycas/force_include.h"
#if defined(POL_PHYCAS)
#	include "pyphy/src/basic_lot.hpp"
#else
#	include "phycas/rand/lot.hpp"
#endif
#include "pyphy/src/probability_distribution.hpp"
#include "pyphy/src/basic_tree.hpp"
#include "pyphy/src/tree_manip.hpp"
//#include "ncl/misc/string_extensions.hpp"
//#include "ncl/nxs_exception.hpp"
#include "pyphy/src/cond_likelihood.hpp"
#include "pyphy/src/cond_likelihood_storage.hpp"
#include "pyphy/src/tip_data.hpp"
#include "pyphy/src/internal_data.hpp"
#include <boost/python.hpp>

#include "pyphy/src/xphylogeny.hpp"

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
		.def("isRoot", &TreeNode::IsRoot)
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
		;

	class_<TreeManip>("TreeManipBase", init<boost::shared_ptr<Tree> >())
		.def("starTree", &TreeManip::starTree)
		.def("randomTree", &TreeManip::randomTree)
		.def("setRandomEdgeLens", &TreeManip::setRandomEdgeLens)
		;

	register_exception_translator<XPhylogeny>(&translateXPhylogeny);
}
