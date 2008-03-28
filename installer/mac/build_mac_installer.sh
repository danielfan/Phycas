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
# python_version=2.3
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
############################  Directory setup  ###################################
#########################################################################

# Specify variables representing the temporary directory wherein packages will be created
stagedir="$PHYCAS_ROOT/installer/mac/stage"
subpackagedir="$stagedir/SubPackages"
boostlibdir="$stagedir/boostlib"
infodir="$PHYCAS_ROOT/installer/mac/info"
packagemakerdir="$PHYCAS_ROOT/installer/mac/packagemaker"
diskimagedir="$PHYCAS_ROOT/installer/mac/diskimage"
userid=$USER
    
# Create a string representing the name of the metapackage
# For example: Phycas-Mac-i386-python2.5.mpkg
buildIdentifier="Phycas-Mac-`arch`-python$python_version"
metapackage="$buildIdentifier.mpkg"

# Recursively delete the "$PHYCAS_ROOT/dist" directory (if it exists). 
rm -rf "$PHYCAS_ROOT/dist" || exit

# Recursively delete the "$PHYCAS_ROOT/build" directory (if it exists). 
rm -rf "$PHYCAS_ROOT/build" || exit

# Recursively delete the "stage" directory (if it exists). 
# This directory is located inside "$PHYCAS_ROOT/phycas/installer/mac"
sudo chown -R $userid.$userid $stagedir
rm -rf "$stagedir" || exit
 
# Now create a new, empty "SubPackages" directory inside the "stage" directory.
# The -p option tells mkdir to create intermediate directories if needed
mkdir -p $subpackagedir || exit
    
# Create a new, empty "boostlib" directory inside the "stage" directory.
# The -p option tells mkdir to create intermediate directories if necessary.
mkdir -p $boostlibdir || exit

echo "********************************************************************************"
echo "Directory setup completed:"
echo "  \$stagedir ($stagedir) deleted and then recreated"
echo "  \$subpackagedir ($subpackagedir) created"
echo "  \$boostlibdir ($boostlibdir) created"
echo "  \$infodir points to $infodir"
echo "  \$packagemakerdir points to $packagemakerdir"
echo "  \$metapackage will be naed $metapackage"
echo "********************************************************************************"
#read -p "Press return to invoke bjam and build phycas..."

#########################################################################
#########################  Invoke bjam  #################################
#########################################################################

# Do the real build of the sources if doBuild equals 1; otherwise, abort.
if test $doBuild = 1
then
    bjam release || exit
fi

echo "********************************************************************************"
echo "Phycas bjam build complete"
echo "********************************************************************************"
#read -p "Press return to invoke packagemaker to build libBoost.pkg..."

#########################################################################
######################  Create boost package  ###########################
#########################################################################

# Copy the boost dylib to the "$boostlib" directory
cp "$PHYCAS_ROOT/phycas/Conversions/$boost_dylib" $boostlibdir

# Now invoke packagemaker to create the boost package:
# -build                         => create an installation package
# -v                             => verbose output during archiving
# -p $subpackagedir/libBoost.pkg       => the path where the package is to be created
# -proj $stagedir/libBoost.pmproj => path to the pmproj document
sudo /Developer/Tools/packagemaker -build -v -p $subpackagedir/libBoost.pkg -proj $packagemakerdir/libBoost.pmproj || exit

echo "********************************************************************************"
echo "Boost dylib copied to $boostlib"
echo "Packagemaker created $subpackagedir/libBoost.pkg"
echo "  using project file $packagemakerlib/libBoost.pmproj"
echo "********************************************************************************"
#read -p "Press return to invoke setup and bdist_mpkg to build Phycas.pkg..."

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

# Copy the current README and LICENSE files to the "$PHYCAS_ROOT/phycas/installer/mac" directory
cp $PHYCAS_ROOT/README "$stagedir/Readme.txt" || exit
cp $PHYCAS_ROOT/LICENSE "$stagedir/License.txt" || exit

# bdist_mpkg already created the phycas package, so just copy it to the SubPackages directory
mv dist/Phycas-*-py$python_version-macosx10*mpkg/Contents/Packages/Phycas-platlib-*-py$python_version-macosx10*.pkg $subpackagedir/Phycas.pkg || exit

echo "********************************************************************************"
echo "distutils.setup used to set up directory structure in $PHYCAS_ROOT/build"
echo "bdist_mpkg used to build Phycas.pkg in $PHYCAS_ROOT/dist"
echo "README and LICENSE files copied to $stagedir as Readme.txt and License.txt, respectively"
echo "Phycas.pkg copied to $subpackagedir"
echo "********************************************************************************"
#read -p "Press return to invoke packagemaker to build the metapackage $metapackage..."

#########################################################################
########################  Create metapackage  ###########################
#########################################################################

# Invoke packagemaker to create a metapackage that will install phycas to site-packages and boost to /usr/local/lib
sudo /Developer/Tools/packagemaker -build -v -p "$stagedir/$metapackage" -v -proj $packagemakerdir/phycas+boost-`arch`-$python_version.pmproj || exit

echo "********************************************************************************"
echo "packagemaker project $stagedir/phycas+boost-`arch`-$python_version.pmproj used to build $stagedir/$metapackage"
echo "********************************************************************************"
#read -p "Press return to create Info.plist..."

#########################################################################
########################  Create Info.plist  ############################
#########################################################################

# Create the correct Info.plist
cp $infodir/info-pList-prefix.txt  $stagedir/info.plist || exit
cat $infodir/10.3test.xml >> $stagedir/info.plist || exit
if test `arch` = "i386" 
then
	cat $infodir/intelTest.xml >> $stagedir/info.plist || exit
else
	cat $infodir/ppcTest.xml >> $stagedir/info.plist || exit
fi
cat $infodir/py${python_version}Test.xml >> $stagedir/info.plist || exit
cat $infodir/info-pList-suffix.txt >> $stagedir/info.plist || exit

# Replace with the correct Info.plist
sudo mv $stagedir/info.plist "$stagedir/$metapackage/Contents/Info.plist" || exit

echo "********************************************************************************"
echo "Info.plist in $stagedir/$metapackage replaced"
echo "********************************************************************************"
#read -p "Press return to create disk image..."

#########################################################################
########################  Create dmg file  ##############################
#########################################################################

# Create the dmg file
cd $diskimagedir || exit
export PHYCAS_MPKG_NAME="$metapackage" || exit
make -f $diskimagedir/dmgMakefile || exit

mv Phycas-1.0.dmg "$PHYCAS_ROOT/dist/${buildIdentifier}.dmg" || exit

echo "********************************************************************************"
echo "$PHYCAS_ROOT/dist/${buildIdentifier}.dmg created"
echo "********************************************************************************"
read -p "Press return to quit..."
