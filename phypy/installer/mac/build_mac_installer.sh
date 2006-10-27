#!/bin/sh
# This script goes to PHYCAS_ROOT and builds a distribution 
#   for your flavor of mac/python
# All temporaries are put on $PHYCAS_ROOT/build and $PHYCAS_ROOT/dist
set -x

if ! test -d "$PHYCAS_ROOT"
then
    echo "PHYCAS_ROOT must be defined"
    exit 1
else
    cd "$PHYCAS_ROOT"
fi
doBuild=1

#do the real build of the source
if test $doBuild = 1
then
    python phypy/dojam.py || exit
fi

# where the packages will go
packPar="$PHYCAS_ROOT/phypy/installer/mac"
subPack="$packPar/SubPackages"
    
mkdir -p "$packPar/lib"
cp "phypy/bin/boost/libs/python/build/libboost_python.dylib/darwin/$BUILD/shared-linkable-true/libboost_python.dylib" "$packPar/lib"

if test -e "$PYTHON_ROOT/bin/bdist_mpkg" ;
then
    #  move the objects and py files to the correct temp location (in build)
    python setup.py build --compiler=fake

    #  build the phycas installer package
    #  tell it that we have compiled.  -k saves temps
    "$PYTHON_ROOT/bin/bdist_mpkg" --skip-build -k
    

    buildIdentifier="Phycas-Mac-`arch`-python$PYTHON_VERSION"
    pack="$buildIdentifier.mpkg"
    cp README "$packPar/Readme.txt"
    cp LICENSE "$packPar/License.txt"
    
    sudo rm -rf "$subPack"
    if ! test -e $subPack
    then
        sudo mkdir -p $subPack || exit
    fi
    
    sudo /Developer/Tools/packagemaker -build -v -p $subPack/libBoost.pkg -proj $packPar/libBoost.pmproj || exit
    
    sudo cp -r dist/Phycas-1.0-py$PYTHON_VERSION-macosx*mpkg/Contents/Packages/Phycas-platlib-1.0-py$PYTHON_VERSION-macosx10*.pkg $subPack/Phycas.pkg || exit
    
    sudo /Developer/Tools/packagemaker -build -v -p "$packPar/$pack" -v -proj $packPar/phycas+boost-`arch`-${PYTHON_VERSION}.pmproj || exit
    
    
    
    cd $packPar || exit

    # create the correct Info.plist
    cp info-pList-prefix.txt  info.plist
    cat 10.3test.xml >> info.plist
    if test `arch` = "i386" 
    then
        cat intelTest.xml >> info.plist
    else
        cat ppcTest.xml >> info.plist    
    fi
    cat py${PYTHON_VERSION}Test.xml >> info.plist
    cat info-pList-suffix.txt >> info.plist

    # replace with the correct Info.plist
    sudo mv info.plist "$pack/Contents/Info.plist"

    echo "$pack created"
    
    #create the dmg file
    cd $packPar
    export PHYCAS_MPKG_NAME="$pack"
    make -f dmgMakefile || exit

    mv Phycas-1.0.dmg "$PHYCAS_ROOT/${buildIdentifier}.dmg"
    
    echo "$PHYCAS_ROOT/${buildIdentifier}.dmg created"
else
    echo "bdist_mpkp is required.  see http://cheeseshop.python.org/pypi/bdist_mpkg/"
    exit 1
fi
