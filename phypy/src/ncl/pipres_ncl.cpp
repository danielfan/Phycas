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

#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/nxs_token.hpp"
#include "ncl/nxs_exception.hpp"
#if defined(PY_PHYCAS)
@
	#include <boost/python/class.hpp>
	#include <boost/python/module.hpp>
	#include <boost/python/def.hpp>
	#include <boost/python/exception_translator.hpp>
	
		// modelled after http://www.boost.org/libs/python/doc/v2/exception_translator.html
	void nxsExceptionToPython(NxsException const & e)
		{ // Use the Python 'C' API to set up an exception object
		PyErr_SetString(PyExc_RuntimeError, e.msg.c_str());
		}	
	
	BOOST_PYTHON_MODULE(pipresNCL)
		{
		using namespace boost::python;

		register_exception_translator<NxsException>(&nxsExceptionToPython);

		class_<NxsToken, bases<>, NxsToken, boost::noncopyable >("NxsToken", no_init)
			// Add a regular member function.
			.def("readNextToken", &NxsToken::ReadToken, return_internal_reference<>())
			;
		
		}

#endif //defined(PHYCAS_EXTENDING_PYTHON)
