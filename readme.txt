15-Oct-2006
-----------
This file is a little out of date now.

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

