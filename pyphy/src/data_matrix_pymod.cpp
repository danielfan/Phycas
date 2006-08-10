#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#include "phycas/force_include.h"
#include "CipresCommlib/CipresDataMatrixHelper.h"

#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(_DataMatrixBase)
{
	class_<CipresNative::DiscreteMatrix, boost::noncopyable>("DataMatrixBase", no_init)
		.def("getNChar", &CipresNative::DiscreteMatrix::getNChar)
		.def("getNTax", &CipresNative::DiscreteMatrix::getNTax)
		.def("getNStates", &CipresNative::DiscreteMatrix::getNStates)
		.def("getSymbolsList", &CipresNative::DiscreteMatrix::getSymbolsList)
		.def("getStateList", &CipresNative::DiscreteMatrix::getStateList, return_value_policy<copy_const_reference>())
		.def("getStateListPos", &CipresNative::DiscreteMatrix::getStateListPos, return_value_policy<copy_const_reference>())
		.def("getCodedDataMatrix", &CipresNative::DiscreteMatrix::getConstNativeC, return_value_policy<copy_const_reference>())
		;
}
