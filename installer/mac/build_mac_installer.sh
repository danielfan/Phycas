#!/bin/sh

# This script goes to $PHYCAS_ROOT and builds a distribution 
# for your flavor of mac/python. Temporaries are put in 
# $PHYCAS_ROOT/build and the final installer ends up in 
# $PHYCAS_ROOT/dist

##########################################################################
#########################  Preliminaries  ################################
##########################################################################

# Use doBuild=1 to rebuild phycas using bjam; use doBuild=0 to skip.
doBuild=1

# Change name of boost dylib below if boost version has changed
boost_dylib="libboost_python-1_34.dylib"

# Set variable python_version to current major.minor version
python_version=`python -c "import platform; major,minor,patch = platform.python_version_tuple(); print '%s.%s' % (major,minor)"`

# Turn tracing on (set -x) or off (set +x). Tracing displays each line 
# as it executes it, showing any substitutions performed.
set -x

##########################################################################
#########################  Sanity Checks  ################################
##########################################################################

# Check to make sure environmental variable $PHYCAS_ROOT is defined; abort if not
if ! test -d "$PHYCAS_ROOT"
then
    echo "PHYCAS_ROOT must be defined"
    exit 1
else
    cd "$PHYCAS_ROOT"
fi

# Check to make sure environmental variable $PYTHON_ROOT is defined; abort if not
if ! test -d "$PYTHON_ROOT"
then
    echo "PYTHON_ROOT must be defined"
    exit 1
fi

# Check if bdist_mpkg exists; if not, abort with message about where to go to get it
if ! test -e "$PYTHON_ROOT/bin/bdist_mpkg" ;
then
    echo "bdist_mpkp is required. See http://cheeseshop.python.org/pypi/bdist_mpkg/"
    exit 1
fi

#########################################################################
#######################  Create temp dirs  ##############################
#########################################################################

# Specify variables representing the temporary directory wherein packages will be created
packPar="$PHYCAS_ROOT/installer/mac"
subPack="$packPar/SubPackages"
    
# Recursively delete the "$PHYCAS_ROOT/dist" directory (if it exists). 
sudo rm -rf "$PHYCAS_ROOT/dist"

# Recursively delete the "$PHYCAS_ROOT/build" directory (if it exists). 
sudo rm -rf "$PHYCAS_ROOT/build"

# Recursively delete the "SubPackages" directory (if it exists). 
# This directory is located inside "$PHYCAS_ROOT/phycas/installer/mac"
sudo rm -rf "$subPack"

# Recursively delete the "lib" directory (if it exists). 
# This directory is located inside "$PHYCAS_ROOT/phycas/installer/mac"
sudo rm -rf "$packPar/lib"

# Now create a new, empty "SubPackages" directory. The -p option tells
# mkdir to create intermediate directories if needed
if ! test -e $subPack
then
	sudo mkdir -p $subPack || exit
fi
    
# Create a directory "lib" inside the "installer/mac" directory. The -p option tells
# mkdir to create intermediate directories if necessary.
mkdir -p "$packPar/lib"

#read -p "temp dirs refreshed"

#########################################################################
#########################  Invoke bjam  #################################
#########################################################################

# Do the real build of the sources if doBuild equals 1; otherwise, abort.
if test $doBuild = 1
then
    bjam release || exit
fi

#########################################################################
######################  Create boost package  ###########################
#########################################################################

# Copy the boost dylib to the "$PHYCAS_ROOT/installer/mac/lib" directory
cp "$PHYCAS_ROOT/phycas/Conversions/$boost_dylib" "$packPar/lib"

# Now invoke packagemaker to create the boost package:
# -build                         => create an installation package
# -v                             => verbose output during archiving
# -p $subPack/libBoost.pkg       => the path where the package is to be created
# -proj $packPar/libBoost.pmproj => path to the pmproj document
sudo /Developer/Tools/packagemaker -build -v -p $subPack/libBoost.pkg -proj $packPar/libBoost.pmproj || exit

#read -p "boost package should be in SubPackages dir now"

#########################################################################
######################  Create phycas package  ##########################
#########################################################################

# Move the objects and py files to the correct temp location (in build).
# --compiler=fake => tells distutils setup() that it should use our
#                    own "compiler" defined in $PHYCAS_ROOT/fakecompiler.py,
#                    which doesn't actually do any compiling but simply 
#                    copies the compiled files to the right location.
# The product of this step is a directory inside $PHYCAS_ROOT/build named
# lib.macosx-10.3-fat-2.5 (the 10.3 comes from a call of the function
# distutils.sysconfig.get_config_var(MACOSX_DEPLOYMENT_TARGET) and 
# represents the minimum version of MacOSX in which this distribution
# could work). Inside lib.macosx-10.3-fat-2.5 is a single directory named
# phycas, which will be the folder that ends up in the user's Python
# site-packages directory.
# 
# The first thing done by setup() is to copy the *.py files in the packages
# listed in the setupArgs dictionary, defined in setup.py. Next, files  
# listed in the data_all list (also defined in setup.py) are copied. Finally,
# fakecompiler.py is invoked, which results in the *.so files for each 
# package being copied into their respective folders in 
# $PHYCAS_ROOT/build/lib.macosx-10.3-fat-2.5/phycas. 
python setup.py build --compiler=fake || exit

# Build the phycas installer package. The --skip-build option tells bdist_mpkg 
# that we have already compiled the sources. The -k option saves temps.
"$PYTHON_ROOT/bin/bdist_mpkg" --skip-build -k || exit

# Create a string representing the name of the metapackage
# For example: Phycas-Mac-i386-python2.5.mpkg
buildIdentifier="Phycas-Mac-`arch`-python$python_version"
pack="$buildIdentifier.mpkg"

# Copy the current README and LICENSE files to the "$PHYCAS_ROOT/phycas/installer/mac" directory
cp $PHYCAS_ROOT/README "$packPar/Readme.txt" || exit
cp $PHYCAS_ROOT/LICENSE "$packPar/License.txt" || exit

# bdist_mpkg already created the phycas package, so just copy it to the SubPackages directory
sudo cp -r dist/Phycas-*-py$python_version-macosx10*mpkg/Contents/Packages/Phycas-platlib-*-py$python_version-macosx10*.pkg $subPack/Phycas.pkg || exit

#read -p "phycas package should be in $subPack dir now"

#########################################################################
########################  Create metapackage  ###########################
#########################################################################

# Invoke packagemaker to create a metapackage that will install phycas to site-packages and boost to /usr/local/lib
sudo /Developer/Tools/packagemaker -build -v -p "$packPar/$pack" -v -proj $packPar/phycas+boost-`arch`-$python_version.pmproj || exit

#read -p "metapackage should be in $packPar dir now"

#########################################################################
########################  Create Info.plist  ############################
#########################################################################

cd $packPar || exit

# Create the correct Info.plist
cp info-pList-prefix.txt  info.plist || exit
cat 10.3test.xml >> info.plist || exit
if test `arch` = "i386" 
then
	cat intelTest.xml >> info.plist || exit
else
	cat ppcTest.xml >> info.plist || exit
fi
cat py${python_version}Test.xml >> info.plist || exit
cat info-pList-suffix.txt >> info.plist || exit

# Replace with the correct Info.plist
sudo mv info.plist "$pack/Contents/Info.plist" || exit

echo "$pack created"

#########################################################################
########################  Create dmg file  ##############################
#########################################################################

# Create the dmg file
cd $packPar || exit
export PHYCAS_MPKG_NAME="$pack" || exit
make -f dmgMakefile || exit

mv Phycas-1.0.dmg "$PHYCAS_ROOT/dist/${buildIdentifier}.dmg" || exit

echo "$PHYCAS_ROOT/dist/${buildIdentifier}.dmg created"

