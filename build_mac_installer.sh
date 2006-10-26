#!/bin/sh
# Steps to creating a  distribution:
#   1. invoke  phypy/dojam.py
#   2. Copy libboost_python.dylib to phypy/phypy/Conversions 
#   3. Make sure "package_data=windows_package_data" is uncommented in setup.py
#   4. Change --target-version below to reflect correct Python target version
#   5. Run this batch file to generate installer (which will be in dist folder)
set -x
cleanOld=1
doBuild=1

if test $cleanOld = 1
then
    sudo rm -rf build
    sudo rm -rf dist
fi

#do the real build of the source
if test $doBuild = 1
then
    python phypy/dojam.py || exit
fi

mkdir -p dist/usr/local/lib
cp "phypy/bin/boost/libs/python/build/libboost_python.dylib/darwin/$BUILD/shared-linkable-true/libboost_python.dylib" dist/usr/local/lib
cp phypy/installer/mac/libBoost.pmproj dist
cp phypy/installer/mac/phycas+boost.pmproj dist

if test -e "$PYTHON_ROOT/bin/bdist_mpkg" ;
then
    python setup.py build --compiler=fake
    "$PYTHON_ROOT/bin/bdist_mpkg" --skip-build -k
    cp README dist/Readme.txt
    cp LICENSE dist/License.txt
    cd dist
    pack=Phycas-Mac-`arch`-python$PYTHON_VERSION.mpkg
    subPack=SubPackages  #$pack/Contents/Packages
    if ! test -e $subPack
    then
        sudo mkdir -p $subPack || exit
    fi
    sudo /Developer/Tools/packagemaker -build -v -p $subPack/libBoost.pkg -proj libBoost.pmproj || exit
    sudo cp -r Phycas-1.0-py$PYTHON_VERSION-macosx*mpkg/Contents/Packages/Phycas-platlib-1.0-py$PYTHON_VERSION-macosx10*.pkg $subPack/ || exit
    sudo /Developer/Tools/packagemaker -build -v -p $pack  -mi Contents -v -proj phycas+boost.pmproj || exit
    echo "$pack created"
else
    echo "bdist_mpkp is required.  see http://cheeseshop.python.org/pypi/bdist_mpkg/"
    exit 1
fi
