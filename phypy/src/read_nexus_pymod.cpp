#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

//#include "phycas/force_include.h"

/*  This file inclusion avoids a bizarre anonymous namespace multiple definition link error that 
	TL (and others) are getting on Mac 10.3.9 (gcc 3.3).
	If you undef PYPHY_INCLUDE_TO_AVOID_LINK_ERROR you need to uncomment:
		cipres_services/cipres_nexus_reader.cpp
		ncl/characters/nxs_charcters_block.cpp
		ncl/command/nxs_command_output.cpp 
	lines from the phypy/read_nexus/Jamfile
*/
#if defined(__APPLE__) && defined(__GNUC__) &&  (__GNUC__ == 3) && (__GNUC_MINOR__ == 3)	//17-May-2006
#	define PYPHY_INCLUDE_TO_AVOID_LINK_ERROR
#endif

#if defined(PYPHY_INCLUDE_TO_AVOID_LINK_ERROR)
#	warning using macro to include cipres_nexus_reader.cpp TEMPORARY HACK!
#	include "phypy/src/cipres/cipres_nexus_reader.cpp"
#	warning using macro to include nxs_charcters_block.cpp TEMPORARY HACK!
#	include "phypy/src/oldphycas/nxs_characters_block.cpp"
#	warning using macro to include nxs_command_output.cpp TEMPORARY HACK!
#	include "phypy/src/ncl/command/nxs_command_output.cpp"
#else
#	include "phypy/src/cipres/cipres_nexus_reader.hpp"
#endif
#include "phypy/src/oldphycas/characters_manager.hpp"
#include "phypy/src/cipres/CipresDataMatrixHelper.h"
#include "phypy/src/ncl/nxs_exception.hpp"
#include "phypy/src/ncl/trees/full_tree_description.hpp"
#include "phypy/src/ncl/misc/nxs_index_set.hpp"

// These next three lines required, otherwise get link error "unresolved external symbol _PyArrayHandle" 
// (at least in VC7)
#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#include <boost/python/numeric.hpp>
#include "phypy/src/thirdparty/num_util/num_util.h"

#include <boost/python.hpp>

using namespace boost::python;


int	CipresNexusReader::GetNChar()
	{
	return phoCharactersMgr->GetNumChars();
	}

void translate(const NxsException & e)
	{
    // Use the Python 'C' API to set up an exception object
    PyErr_SetString(PyExc_Exception, e.what());
    }

BOOST_PYTHON_MODULE(_ReadNexus)
{
	class_<FullTreeDescription>("FullTreeDescription", init<FullTreeDescription>())
		.def_readwrite("name", &FullTreeDescription::name)
		.def_readwrite("newick", &FullTreeDescription::newick)
		.def_readwrite("rooted", &FullTreeDescription::rooted)
		;

#if 1
	// If this section is included, Visual Studio .NET 2003 Service Pack 1 produces
	// error C2247: 'boost::python::objects::iterator_range<NextPolicies,Iterator>::next' not accessible because 'boost::details::compressed_pair_imp<T1,T2,Version>' uses 'private' to inherit from 'boost::python::objects::iterator_range<NextPolicies,Iterator>::next'
	class_<std::vector<FullTreeDescription> >("VecFullTreeDescription", no_init)
		.def("__iter__",  iterator<std::vector<FullTreeDescription> >())
		;
#endif
   
	class_<CipresNexusReader>("NexusReaderBase", init<int>())
        .def("readFile", &CipresNexusReader::ReadFilePath)
        .def("getNChar", &CipresNexusReader::GetNChar)
        .def("getErrorMessage", &CipresNexusReader::GetErrorMessage)
        .def("getTrees", &CipresNexusReader::GetTrees, return_value_policy<copy_const_reference>())
		.def("getTaxLabels", &CipresNexusReader::GetTaxLabels)
		;
	
	def("getDiscreteMatrix", createNativeDiscreteMatrix, return_value_policy<manage_new_object>());

	register_exception_translator<NxsException>(&translate);
}
