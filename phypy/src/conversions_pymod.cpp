#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle

//#include "phycas/force_include.h"
#include "phypy/src/basic_tree.hpp"
#include <boost/python.hpp>
//#include <boost/python/tuple.hpp>
//#include <boost/python/numeric.hpp>
#include "phypy/src/thirdparty/num_util.h"
#include "phypy/src/thirdparty/pyconversions.h"	// from HippoDraw
#include "phypy/src/cipres/ConfigDependentHeaders.h"	// int8_t typedef //POL 18Mar2006
using namespace boost::python;

BOOST_PYTHON_MODULE(_Conversions)
{
	// these lines required by num_util
	import_array();
	numeric::array::set_module_and_type(); // defaults to "numarray", "NDArray"
	//numeric::array::set_module_and_type("Numeric", "ArrayType");	// old numeric, don't use

	// these lines taken from HippoGraph
	std_vector_to_tuple<unsigned>();
	std_vector_to_tuple<int8_t>();
	std_vector_to_tuple<double>();
	std_vector_to_tuple<std::string>();
	std_vector_to_tuple<phycas::TreeNode *>();

	from_python_sequence<std::vector<unsigned>, variable_capacity_policy>();
    from_python_sequence<std::vector<int8_t>, variable_capacity_policy>();
    from_python_sequence<std::vector<double>, variable_capacity_policy>();
    from_python_sequence<std::vector<std::string>, variable_capacity_policy>();
}
