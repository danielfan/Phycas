First, note that the *proper* way to build Phycas is using bjam from the directory containing
the file Jamroot (i.e. the parent of this directory). There is really only one good reason 
you should be trying to use this Visual Studio solution: you are interested in debugging C++
Phycas code using the IDE. If this is the case, you will need to do some tweaking before this
"solution" will work:

Setting environmental variables
-------------------------------
You will need to set some environmental variables before building Phycas using this 
Visual Studio solution:

BOOST_ROOT=C:\boost_1_34_0

Building debug Python 2.5 library
---------------------------------
You need python25_d.lib in order to build a debug verion of Phycas, but unfortunately
python25_d.lib does not come with the standard Python 2.5 release - you have to download
the Python source code and build it yourself. If you try to build a debug version of 
Phycas without having this debug version of Python, you will get link errors such as 
these:

phylogeny_pymod.obj : error LNK2019: unresolved external symbol __imp___Py_Dealloc ...
phylogeny_pymod.obj : error LNK2019: unresolved external symbol __imp___Py_NegativeRefcount ...
phylogeny_pymod.obj : error LNK2001: unresolved external symbol __imp___Py_RefTotal ...


Building the Boost Python debug DLL
-----------------------------------
To build *only* the Python library out of the many projects composing Boost, navigate
to the Boost folder (e.g. C:\boost_1_34_0) and type

bjam debug --with-python stage
             
This will cause the boost python dll and link lib to be placed in C:\boost_1_34_0\stage\lib.
                                              
C:\Python25\include\pyconfig.h
------------------------------
Python, in trying to be helpful, nominates python25_d.lib as the Python link library to be used
when in building a debug verison. Unfortunately, python25_d.lib is not distributed with Python,
so this isn't all that helpful after all. Here is the offending code (starting on line 275 in 
Python25\include\pyconfig.h):

	/* For an MSVC DLL, we can nominate the .lib files used by extensions */
	#ifdef MS_COREDLL
	#	ifndef Py_BUILD_CORE /* not building the core - must be an ext */
	#		if defined(_MSC_VER)
				/* So MSVC users need not specify the .lib file in
				their Makefile (other compilers are generally
				taken care of by distutils.) */
	#			ifdef _DEBUG
	#				pragma comment(lib,"python25_d.lib")
	#			else
	#				pragma comment(lib,"python25.lib")
	#			endif /* _DEBUG */
	#		endif /* _MSC_VER */
	#	endif /* Py_BUILD_CORE */
	#endif /* MS_COREDLL */

The rather inelegant way that I deal with this is to just disable the section beginning with 
"#if defined(_MSC_VER)":

	/* For an MSVC DLL, we can nominate the .lib files used by extensions */
	#ifdef MS_COREDLL
	#	ifndef Py_BUILD_CORE /* not building the core - must be an ext */
	#		if 0 && defined(_MSC_VER)	//POL disabled
				/* So MSVC users need not specify the .lib file in
				their Makefile (other compilers are generally
				taken care of by distutils.) */
	#			ifdef _DEBUG
	#				pragma comment(lib,"python25_d.lib")
	#			else
	#				pragma comment(lib,"python25.lib")
	#			endif /* _DEBUG */
	#		endif /* _MSC_VER */
	#	endif /* Py_BUILD_CORE */
	#endif /* MS_COREDLL */

C:\boost_1_34_0\boost\python\detail\caller.hpp
----------------------------------------------
Cast the return value on line 52 of this file to avoid many warnings while compiling Phycas.
Before:
    return PyTuple_GET_SIZE(args_);
After:
    return (unsigned)PyTuple_GET_SIZE(args_);

C:\boost_1_34_0\boost\python\detail\config.hpp
----------------------------------------------
At the end of this file, Boost Python also tries to be helpful:

	//  enable automatic library variant selection  ------------------------------// 
	
	#if !defined(BOOST_PYTHON_SOURCE) && !defined(BOOST_ALL_NO_LIB) && !defined(BOOST_PYTHON_NO_LIB)
	//
	// Set the name of our library, this will get undef'ed by auto_link.hpp
	// once it's done with it:
	//
	#define BOOST_LIB_NAME boost_python
	//
	// If we're importing code from a dll, then tell auto_link.hpp about it:
	//
	#ifdef BOOST_PYTHON_DYNAMIC_LIB
	#  define BOOST_DYN_LINK
	#endif
	//
	// And include the header that does the work:
	//
	#include <boost/config/auto_link.hpp>
	#endif  // auto-linking disabled

As a result, you get this link error when trying to compile Phycas:

	LINK : fatal error LNK1104: cannot open file 'boost_python-vc71-mt-gd-1_34.lib'

The library name is boost_python_debug.lib when the VC project supplied in 
C:\boost_1_34_0\libs\python\build\VisualStudio is used to compile the Boost Python DLLs.
I have defined BOOST_PYTHON_SOURCE in the Phycas VC projects to avoid this.
