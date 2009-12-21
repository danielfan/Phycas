Last updated 11 Dec 2009 by Paul Lewis

XCode project for Phycas
------------------------

An XCode project for Phycas was created in the first week of December, 2009, as an alternative
to the bjam build system. Having an XCode project is especially useful for compiling because
bjam takes some time to determine what needs to be recompiled, whereas in XCode one can explicitly
recompile just one file if desired, which is much faster. 

What this XCode project does (and does not do) for you
------------------------------------------------------
This XCode project builds libboost_python.dylib and libncl.dylib, as well as the following files:
_ConversionsExt.so, _DataMatrixExt.so, _LikelihoodExt.so, _PhylogenyExt.so, _ProbDistExt.so, and
_ReadNexusExt.so. It does not build a debug version of Python, so you must do this by downloading
the Python tgz file, uncompressing it, then running configure and make to obtain libpython2.6.dylib.
Here is the sequence I used:

cd /Users/plewis/Python-2.6.4
./configure --with-pydebug --enable-shared --prefix=/Users/plewis/pydbg
make install

Afterwards, the libpython2.6.dylib file will be located in /Users/plewis/pydbg/lib.

Setting up a custom executable
------------------------------
To run a phycas Python script, you need to first set up a custom executable within XCode. These
custom executables need to be defined by every user separately (hence "custom"). Right-click on the 
"Executables" part of the project and choose "Add|New Custom Executable..." Choose a name for your
custom executable (e.g. "run_python), replacing the default name ("Executable"), and use the browse 
button to select the debug version of the python executable (e.g. /Users/plewis/pydbg/bin/python if 
you used my instructions above to build Python). After clicking the Next button to add this to your
phycas project, you will get a 4-part dialog that has tabs named "General", "Arguments", "Debugging"
and "Comments". In the "General" pane, set the working directory to "Custom directory" and specify
the absolute path to the directory in which the script you wish to run resides. In the "Arguments"
pane, specify the name of the script you want to run as an argument, and specify the following 
environmental variables in the lower part:
	DYLD_LIBRARY_PATH	$(PHYCAS_ROOT)/phycas/Conversions
	PYTHONPATH			$(PHYCAS_ROOT)
Finally, in the "Debugging" section, be sure "Break on Debugger() and DebugStr()" is checked.

Lazy symbol loading
-------------------
Because of the somewhat non-standard nature of phycas (combining C++ with Python code), it is 
necessary to uncheck "Load symbols lazily" in the "symbol Loading Options" section of the
Debugging tab in the main XCode Preferences dialog. Failing to do with will result in breakpoints
that are ignored (although DebugStr calls still work: see below) and turn yellow when the program
is running.

Using DebugStr
--------------
Although using ordinary breakpoints are usually preferable, you can include code to cause the 
debugger to stop at a place of your choosing. (Note that you might not realize that the debugger 
has stopped, for it doesn't pop to the foreground automatically.) For example, including this line...

  DebugStr("\pin harvestLnLFromValidNode");
  
will stop the debugger after displaying the string "in harvestLnLFromValidNode" (the initial \p
is required because this is a "Pascal string"). At the top of a file in which you want to 
use DebugStr, include the following

#include <CoreServices/CoreServices.h>
#undef check	

The "#undef check" is needed because CoreServices.h defines a macro named "check" and this 
conflicts with some phycas functions (e.g. UnderflowManager::check).

Drawbacks
---------
Unfortunately, there seems to be one bug in XCode that makes sharing it across different users
a bit inconvenient at present. Hopefully we can find a workaround in the near future for this 
problem.

Even though we specified Path Type to be "Relative to PYTHON_ROOT" for libpython2.6.dylib, this
one item in the "External Frameworks and Libraries" section of the project needs to be added by
every user separately. If you simply update to obtain a fresh version of the XCode project, you
will find link errors when you build. These will go away if you delete libpython 2.6.dylib from
the project, then simply add it back in by dragging it from a Finder window to the "External 
Frameworks and Libraries" section of the project. We think this is because the file 
project.pbxproj (within phycas.xcodeproj) stores a file reference that is unique to a particular
user (and perhaps a particular computer) and invalid when the user changes. Thus, whenever any
one of us (svn) commits, the XCode project will be screwed up for everyone else when they (svn) 
update.
   
