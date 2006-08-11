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

The following clean out the so/pyd dynamic libraries. This seems to be necessary
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

PYTHON_ROOT     C:\Python24
PYTHON_VERSION  2.4
BOOST_ROOT      C:\boost_1_33_0
PHYCAS_ROOT     C:\Synchronized\Projects\phycasdev

If the bjam build is successful, you will need to add something like this to 
your PYTHONPATH environmental variable:

PYTHONPATH      $HOME/phycasdev/pyphy

To try running the examples (in total, this should take only a couple of minutes):

cd phycasdev/pyphy/pyphy/Examples
./runall.sh

If the script does not abort, then the output of each example matched the reference
output and all is well.

