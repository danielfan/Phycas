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
