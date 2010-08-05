Last updated 19 Dec 2008 by Paul O. Lewis

First, note that the *proper* way to build Phycas is using bjam from the directory containing
the file Jamroot (i.e. the parent of this directory). There is really only one good reason 
you should be trying to use this Visual Studio solution: you are interested in debugging C++
Phycas code using the IDE. If this is the case, you will need to do some tweaking before this
"solution" will work. 

Setting environmental variables
-------------------------------
You will need to set some environmental variables before building Phycas using this 
Visual Studio solution:

BOOST_ROOT=C:\Users\Administrator\Library\boost_1_37_0
NCL_ROOT=C:\Users\Administrator\Documents\ndev\branches\v2.1
PHYCAS_NCL_STATIC=1
PHYCAS_ROOT=C:\Users\Administrator\Documents\pdev\trunk
PYTHON_ROOT=C:\Python26
PYTHONPATH=%PHYCAS_ROOT%
PATH=%PATH%;C:\Python26;C:\Python-2.6.1\PCbuild;C:\Users\Administrator\Library\boost_1_37_0\tools\jam\src\bin.ntx86

Building debug Python 2.6 library
---------------------------------
You need python_d.exe, python26_d.dll and python26_d.lib in order to build a debug version of Phycas, 
but unfortunately none of these come with the standard Python 2.6 release - you have to download
the Python source code and build it yourself. If you try to build a debug version of 
Phycas without having this debug version of Python, you will get link errors such as 
these:

phylogeny_pymod.obj : error LNK2019: unresolved external symbol __imp___Py_Dealloc ...
phylogeny_pymod.obj : error LNK2019: unresolved external symbol __imp___Py_NegativeRefcount ...
phylogeny_pymod.obj : error LNK2001: unresolved external symbol __imp___Py_RefTotal ...

To build it yourself, go to the directory C:\Python-2.6.1\PCbuild and open the pcbuild.sln file, 
choose Debug version and then build just the "python" project (this is the only project you need to build).
Be sure to add C:\Python-2.6.1\PCbuild to the PATH environmental variable so that python_d.exe can
be found when it is invoked.

Building Tcl/Tk into Python debug build
---------------------------------------
To use Tkinter in debug mode, tried to follow the directions for building the _tkinter subproject in 
C:\Python-2.6.1\PCbuild\readme.txt, but a lot seems to have changed. Here is a synopsis of all I had
to do to compile the _tkinter subproject in C:\Python-2.6.1\PCbuild\pcbuild.sln:

Download tcl856-src.zip and tk856-src.zip from http://www.tcl.tk/software/tcltk/download.html
Extract files from both zip archives to C:\tcl8.5.6 and C:\tk8.5.6
Run C:\tcl8.5.6\win\buildall.vc.bat debug
    I had to edit buildall.vc.bat first to specify correct location for vcvars32.bat
    I got 1 error (could not find hcrtf) that appears to be related to creating help/documentation
Run C:\tk8.5.6\win\buildall.vc.bat debug
    I had to edit buildall.vc.bat first to specify correct location for vcvars32.bat
    Before running, specify "set TCLDIR=..\..\tcl8.5.6"
    I got 1 error (could not find hcrtf) that appears to be related to creating help/documentation
Specify the following for additional include directories in the _tkinter project:
    C:\tcl8.5.6\Generic;C:\tk8.5.6\Generic;C:\tk8.5.6\xlib
Specify the following for additional library directories in the _tkinter project:
    C:\tcl8.5.6\win\Debug_VC9;C:\tk8.5.6\win\Debug_VC9;
Specify the following for additional dependencies in the linker section of the _tkinter project:
    tcl85g.lib tk85g.lib
Build the python project again to regenerate python_d.exe

Compile the _tkinter project in Python-2.6.1\PCBuild folder\sln

Building the Debug version of the Boost Python DLL
--------------------------------------------------
Modify C:\boost_1_37_0\tools\build\v2\user-config.jam to include these lines:

using python : 2.6 : C:\\Python26\\python.exe ;
using python : 2.6 : C:\\Python-2.6.1\\PCBuild\\python_d.exe
  : # includes
  : # libs
  : <python-debugging>on ;

Failure to add these lines results in a boost_python*.dll that loads non-debugging DLLs,
causing python_d.exe to crash. Also, in user-config.jam, make sure the following line is 
uncommented:

using msvc ;

To build *only* the Python library out of the many projects composing Boost, navigate
to the Boost folder (e.g. C:\boost_1_37_0) and type

bjam debug --with-python python-debugging=on stage

This will cause the boost python dll (boost_python-vc90-mt-gyd-1_37.dll) and import lib 
(boost_python-vc90-mt-gyd-1_37.lib) to be placed in C:\boost_1_37_0\stage\lib.

Note that the python-debugging feature is important (although, in theory, it is not
necessary to include "python-debugging=on" in the bjam command line because we have
specified this option in "user-config.jam"). Without python-debugging set to "on", the 
libraries boost generates will be named boost_python-vc90-mt-gd-1_37.lib and 
boost_python-vc90-mt-gd-1_37.dll, and you will encounter an error like this when 
attempting to import phycas within the python_d interpreter:

Python 2.5.1 (r251:54863, May 31 2007, 11:35:25) [MSC v.1310 32 bit (Intel)] on win32
Type "help", "copyright", "credits" or "license" for more information.
>>> import phycas
Fatal Python error: Interpreter not initialized (version mismatch?)

(Obviously, that snippet was from an older Python 2.5 distribution, and I haven't bothered
to update it.) An even worse fate awaits you should you try using the regular release-mode 
python interpreter:

Python 2.5.1 (r251:54863, Apr 18 2007, 08:51:08) [MSC v.1310 32 bit (Intel)] on win32
Type "help", "copyright", "credits" or "license" for more information.
>>> import phycas
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "C:\Synchronized\Projects\phycasdev_trunk\phycas\__init__.py", line 15, in <module>
    import Likelihood
  File "C:\Synchronized\Projects\phycasdev_trunk\phycas\Likelihood\__init__.py", line 5, in <module>

    from _LikelihoodBase import *
RuntimeError: unidentifiable C++ exception
                                              
C:\Python26\include\pyconfig.h
------------------------------
Python, in trying to be helpful, nominates python26_d.lib as the Python link library to be used
when in building a debug verison. Here is the offending code (starting on line 318 in 
Python26\include\pyconfig.h):

	/* For an MSVC DLL, we can nominate the .lib files used by extensions */
	#ifdef MS_COREDLL
	#	ifndef Py_BUILD_CORE /* not building the core - must be an ext */
	#		if defined(_MSC_VER)
				/* So MSVC users need not specify the .lib file in
				their Makefile (other compilers are generally
				taken care of by distutils.) */
	#			ifdef _DEBUG
	#				pragma comment(lib,"python26_d.lib")
	#			else
	#				pragma comment(lib,"python26.lib")
	#			endif /* _DEBUG */
	#		endif /* _MSC_VER */
	#	endif /* Py_BUILD_CORE */
	#endif /* MS_COREDLL */

In the past, I've simply commented out the part where it nominates python26_d.lib, but
that doesn't seem to work anymore. Hence the need to actually build python26_d.lib and
python26_d.dll from source. Note that you do not need to actually do anything to pyconfig.h:
I included this section merely to remind myself why VC knows it should link in python26_d.lib 
when I haven't specified this anywhere in the project!

C:\boost_1_37_0\boost\python\detail\caller.hpp
----------------------------------------------
Cast the return value on line 55 of this file to avoid many warnings while compiling Phycas.
Before:
    return PyTuple_GET_SIZE(args_);
After:
    return (unsigned)PyTuple_GET_SIZE(args_);

C:\boost_1_34_0\boost\python\detail\config.hpp
----------------------------------------------
Nothing to be done here. Just a reminder that the end of this file seems to be where Boost
comes up with the rather complicated name of the import lib to link (e.g. )

Other tweaks I've made in the various projects within Vc.sln
------------------------------------------------------------
C/C++: Additional include directories:

	$(PYTHON_ROOT)\include;$(PHYCAS_ROOT);$(BOOST_ROOT)
	
C/C++: Code generation:
	
	Multi-threaded Debug DLL (/MDd)   [debug version]
	Multi-threaded DLL (/MD)          [release version]
	
C/C++: Language: 

	Enable Run-Time Type Info: Yes (/GR)
	
C/C++: Advanced: Force Includes: 

	$(PHYCAS_ROOT)\phycas\src\std_force_include.hpp
	
Linker: General: Output File: 

	$(PHYCAS_ROOT)/phycas/Conversions/_Conversions.pyd    [release version]
	$(PHYCAS_ROOT)/phycas/Conversions/_Conversions_d.pyd  [debug version]
	
	Note: the "_d" suffix is important. Python will generate errors like the following if you use 
	python_d.exe as your interpreter but fail to add the _d to your extension DLL:
	
		ImportError: No module named _Conversions
	
Linker: General: Additional Library Directories: 

	$(BOOST_ROOT)\stage\lib;$(PYTHON_ROOT)\libs

Linker: Input: Additional Dependencies: 

	boost_python-vc90-mt-gyd-1_37.lib  [variant created by python-debugging property]

Linker: Input: Ignore Specific Library: 

	boost_python-vc90-mt-gd-1_37.lib   [this is the one Boost suggests, but this one will not work]

In addition, anytime C files are included in the project, you need to select their names, choose Properties from
the right-click menu, and set C/C++: Advanced: Force Includes to $(NOINHERIT)

Compiling the extensions
------------------------
In the Solution Explorer, right-click "Solution 'VC' (6 projects)" and choose Build. All six python extensions should
build without errors or warnings with the exception of the third-party DCFLIB file dcdflib.c, which issues the following
warnings:

c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(464) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(479) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(491) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(717) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(922) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(1073) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(1213) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(1303) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(1304) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(1324) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(1325) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(5256) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(6605) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(6620) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(6632) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(6723) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(7896) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data
c:\Synchronized\Projects\phycasdev_trunk\phycas\src\thirdparty\dcdflib\src\dcdflib.c(8080) : warning C4244: '=' : conversion from 'double' to 'int', possible loss of data

I have chosen to not fix these warnings out of a desire to leave this library completely unmodified. 
 
After compilation
-----------------
After the phycas python extensions have compiled, you should be able to start the debug version of python (python_d.exe)
and import phycas. 

Deja vu
-------
o "Compiler out of keys" errors in VC debug builds arise from using __FILE__ and __LINE__ macros with the 
compiler option /FI (). I have changed /FI to /Fi in every project and now the error no longer arises.
Note that even if you do not use __FILE__ and __LINE__ explicitly in your code, these macros are still
being used if you use other macros that employ them (e.g. assert).



