This is the phycas readme file for developers with SVN access who are 
building with Boost's bjam tool (this is the recommended build system for
phycas).

29-May-2007
-----------

Build Prerequisites
-------------------

   * You must have the Boost Build v2 version of bjam compiled
   * BOOST_ROOT must be in your environment and must refer to the top of your
     boost libraries
   * The boost build system must be configured. For instance the
     file  $BOOST_ROOT/tools/build/v2/user-config.jam must contain the correct
     settings
   * You should have one of the following settings to indicate your Operating 
     System:
     	For linux:     OSTYPE=linux
		For mac:       OSTYPE=darwin
        For Windows:   OS=Windows_NT

Building
--------

Simply type the command:
   bjam
   

Prerequisites for Running
-------------------------
If you do not install the phycas directory in a standard location (such as
the site-packages directory in your python library) then you will have to
add the PHYCAS_ROOT directory (the top level directory that contains the "phycas"
directory) to your PYTHONPATH environmental setting.

Testing
-------

Tests can be invoke as shown below (commands after the $ prompt output on the 
next line:

$ cd phycas/Tests
$ python doctestall.py
Finished testing examples
$ python runall.py

The runall.py script will produce a lot of output.

Miscellaneous information
-------------------------

The following file contains the "release" history (what features were added when) 
and shows projections for future feature additions:

  CHANGES

The following files are used by bjam:

  Jamroot




23-Oct-2006
-----------
This file is a little out of date now.

The binary distribution is all under the folder below:

  phypy/phypy 

Each of the following folders contain the Jamfile and VC projects needed 
for building the corresponding module directory of the same name inside 
the "<phycasdev>/pyphy/pyphy" folder:

  phypy/conversions
  phypy/data_matrix
  phypy/likelihood
  phypy/phylogeny
  phypy/prob_dist
  phypy/read_nexus

All C++ source code (and header files) are in:

  phypy/src

The master Visual Studio solution is in the folder below (but not yet working):

  phypy/VC

The following file contains the "release" history (what features were added when) 
and shows projections for future feature additions:

  CHANGES

The following cleans out the so/pyd dynamic libraries. This seems to be necessary
(on Windows at least) when switching between debug and release versions:

  phypy/clean.sh     (shell script)
  phypy/clean.bat   (Windows batch file)

The following initiate the build process using bjam:

  phypy/dojam.sh       (shell script)
  phypy/dojam.bat   (Windows batch file)

The following files are used by bjam:

  Jamfile
  Jamrules

There are several environmental variables that dojam and/or dojam.bat expect to 
be set before they are run. Here are the environmental variable settings used 
on a Windows machine in case it helps:

PYTHON_ROOT     C:\Python25
PYTHON_VERSION  2.5
BOOST_ROOT      C:\boost_1_33_1
PHYCAS_ROOT     C:\Synchronized\Projects\phycasdev

If the bjam build is successful, you will need to add something like this to 
your PYTHONPATH environmental variable:

PYTHONPATH      $HOME/phycasdev/pyphy 

You also need to put boost_python.dll where it can be found. I usually copy it from
  phycasdev_trunk\phypy\bin\boost\libs\python\build\boost_python.dll\vc-7_1\release\threading-multi
to
  phycasdev_trunk\phypy\phypy\Conversions
When you issue the "from phypy import *" command in Python, the first module that is
loaded is Conversions, and that's also the first time the system tries to find
boost_python.dll. If boost_python.dll is in the same directory as _Conversions.pyd,
then the system will find boost_python.dll there too and all will be well

To try running the tests (in total, this should take only a couple of minutes):

cd phycasdev/pyphy/pyphy/Tests
runall.bat

If the script does not abort, then the output of each test run matched the reference
output and all is well.

/src

The master Visual Studio solution is in the folder below (but not yet working):
  phycasdev/pyphy/VC

The following file contains the "release" history (what features were added when) 
and shows projections for future feature additions:

  phycasdev/pyphy/changes.txt

The following cleans out the so/pyd dynamic libraries. This seems to be necessary
(on Windows at least) when switching between debug and release versions:

  phycasdev/pyphy/clean       (shell script)
  phycasdev/pyphy/clean.bat   (Windows batch file)

The following initiate the build process using bjam:

  phycasdev/pyphy/dojam       (shell script)
  phycasdev/pyphy/dojam.bat   (Windows batch file)

The following files are used by bjam:

  Jamfile
  Jamrules

There are several environmental variables that dojam and/or dojam.bat expect to 
be set before they are run. Here are the environmental variable settings used 
on a Windows machine in case it helps:

PYTHON_ROOT     C:\Python25
PYTHON_VERSION  2.5
BOOST_ROOT      C:\boost_1_33_1
PHYCAS_ROOT     C:\Synchronized\Projects\phycasdev

If the bjam build is successful, you will need to add something like this to 
your PYTHONPATH environmental variable:

PYTHONPATH      $HOME/phycasdev/pyphy 

You also need to put boost_python.dll where it can be found. I usually copy it from
  phycasdev_trunk\phypy\bin\boost\libs\python\build\boost_python.dll\vc-7_1\release\threading-multi
to
  phycasdev_trunk\phypy\phypy\Conversions
When you issue the "from phypy import *" command in Python, the first module that is
loaded is Conversions, and that's also the first time the system tries to find
boost_python.dll. If boost_python.dll is in the same directory as _Conversions.pyd,
then the system will find boost_python.dll there too and all will be well

To try running the tests (in total, this should take only a couple of minutes):

cd phycasdev/pyphy/pyphy/Tests
runall.bat

If the script does not abort, then the output of each test run matched the reference
output and all is well.



11-August-2006
--------------
In this directory, you should find only two folders: "pyphy" and "unused". 
The "unused" folder is where everything that isn't currently being used went, 
so if you are looking for something, poke around in there for it. Inside the 
"pyphy" folder you will find these folders and files (the folder that you 
checked out into via SVN is called phycasdev below):

The binary distribution is all under the folder below:

  phycasdev/pyphy/pyphy 

Each of the following folders contain the Jamfile and VC projects needed 
for building the corresponding module directory of the same name inside 
the "<phycasdev>/pyphy/pyphy" folder:

  phycasdev/pyphy/conversions
  phycasdev/pyphy/data_matrix
  phycasdev/pyphy/likelihood
  phycasdev/pyphy/phylogeny
  phycasdev/pyphy/prob_dist
  phycasdev/pyphy/read_nexus

All source code and header files are now under this one roof:

  phycasdev/pyphy/src

The master Visual Studio solution is in the folder below (but not yet working):
  phycasdev/pyphy/VC

The following file contains the "release" history (what features were added when) 
and shows projections for future feature additions:

  phycasdev/pyphy/changes.txt

The following cleans out the so/pyd dynamic libraries. This seems to be necessary
(on Windows at least) when switching between debug and release versions:

  phycasdev/pyphy/clean       (shell script)
  phycasdev/pyphy/clean.bat   (Windows batch file)

The following initiate the build process using bjam:

  phycasdev/pyphy/dojam       (shell script)
  phycasdev/pyphy/dojam.bat   (Windows batch file)

The following files are used by bjam:

  Jamfile
  Jamrules

There are several environmental variables that dojam and/or dojam.bat expect to 
be set before they are run. Here are the environmental variable settings used 
on a Windows machine in case it helps:

PYTHON_ROOT     C:\Python25
PYTHON_VERSION  2.5
BOOST_ROOT      C:\boost_1_33_1
PHYCAS_ROOT     C:\Synchronized\Projects\phycasdev

If the bjam build is successful, you will need to add something like this to 
your PYTHONPATH environmental variable:

PYTHONPATH      $HOME/phycasdev/pyphy 

You also need to put boost_python.dll where it can be found. I usually copy it from
  phycasdev_trunk\phypy\bin\boost\libs\python\build\boost_python.dll\vc-7_1\release\threading-multi
to
  phycasdev_trunk\phypy\phypy\Conversions
When you issue the "from phypy import *" command in Python, the first module that is
loaded is Conversions, and that's also the first time the system tries to find
boost_python.dll. If boost_python.dll is in the same directory as _Conversions.pyd,
then the system will find boost_python.dll there too and all will be well

To try running the tests (in total, this should take only a couple of minutes):

cd phycasdev/pyphy/pyphy/Tests
runall.bat

If the script does not abort, then the output of each test run matched the reference
output and all is well.

