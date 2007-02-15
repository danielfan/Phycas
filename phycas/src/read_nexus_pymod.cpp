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

//#include "phycas/force_include.h"

/*  This file inclusion avoids a bizarre anonymous namespace multiple definition link error that 
	TL (and others) are getting on Mac 10.3.9 (gcc 3.3).
	If you undef PYPHY_INCLUDE_TO_AVOID_LINK_ERROR you need to uncomment:
		cipres_services/cipres_nexus_reader.cpp
		ncl/characters/nxs_charcters_block.cpp
		ncl/command/nxs_command_output.cpp 
	lines from the phycas/read_nexus/Jamfile
*/
#if defined(__APPLE__) && defined(__GNUC__) &&  (__GNUC__ == 3) && (__GNUC_MINOR__ == 3)	//17-May-2006
#	define PYPHY_INCLUDE_TO_AVOID_LINK_ERROR
#endif

#if defined(PYPHY_INCLUDE_TO_AVOID_LINK_ERROR)
#	warning using macro to include cipres_nexus_reader.cpp TEMPORARY HACK!
#	include "phycas/src/cipres/cipres_nexus_reader.cpp"
#	warning using macro to include nxs_charcters_block.cpp TEMPORARY HACK!
#	include "phycas/src/oldphycas/nxs_characters_block.cpp"
#	warning using macro to include nxs_command_output.cpp TEMPORARY HACK!
#	include "phycas/src/ncl/command/nxs_command_output.cpp"
#else
#	include "phycas/src/cipres/cipres_nexus_reader.hpp"
#endif
#include "phycas/src/oldphycas/characters_manager.hpp"
#include "phycas/src/cipres/CipresDataMatrixHelper.h"
#include "phycas/src/ncl/nxs_exception.hpp"
#include "phycas/src/ncl/trees/full_tree_description.hpp"
#include "phycas/src/ncl/misc/nxs_index_set.hpp"

#if defined(USING_NUMARRAY)
// These next three lines required, otherwise get link error "unresolved external symbol _PyArrayHandle" 
// (at least in VC7)
#	define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#	include <boost/python/numeric.hpp>
#	include "phycas/src/thirdparty/num_util/num_util.h"
#endif

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

	class_<std::vector<FullTreeDescription> >("VecFullTreeDescription", no_init)
		.def("__iter__",  iterator<std::vector<FullTreeDescription> >())
		;
   
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
