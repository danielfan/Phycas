export PATH=/usr/local/bin:/usr/local/sbin:/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/teTeX/bin/i386-apple-darwin-current
unset DYLD_LIBRARY_PATH
unset PYTHONPATH

export SELF_CONTAINED_BUILD_ROOT=/Users/mholder/Desktop/stepping_stone_lab/mac/build

if test -z $SELF_CONTAINED_CONTENTS_DIR
then
    if test -z ${SELF_CONTAINED_BUILD_ROOT}
    then
        echo "SELF_CONTAINED_BUILD_ROOT must be defined"
        exit
    fi
    export SELF_CONTAINED_CONTENTS_DIR="${SELF_CONTAINED_BUILD_ROOT}/PhycasGUI.app/Contents"
    export SELF_CONTAINED_PREFIX="${SELF_CONTAINED_CONTENTS_DIR}/Resources"
fi


export DYLD_LIBRARY_PATH="$SELF_CONTAINED_PREFIX/lib:$SELF_CONTAINED_PREFIX/lib/ncl:$SELF_CONTAINED_PREFIX/lib/python2.6/site-packages/phycas:$DYLD_LIBRARY_PATH:"
export CXXFLAGS="-arch i386 -m32"
export CFLAGS="-arch i386 -m32"
export LDFLAGS="-arch i386 -m32"
export MACOSX_DEPLOYMENT_TARGET="10.4"
export BOOST_ROOT="${SELF_CONTAINED_BUILD_ROOT}/boost_1_43_0"
export BOOST_BUILD_PATH="${BOOST_ROOT}/tools/build/v2"
export PHYCAS_ROOT="${SELF_CONTAINED_BUILD_ROOT}/phycas_trunk"
export OSTYPE="darwin"
export PATH="${SELF_CONTAINED_PREFIX}/bin:${BOOST_ROOT}/tools/jam/src/bin.macosxx86_64:$SELF_CONTAINED_PREFIX/bin:$PATH"
export NCL_INSTALL_DIR="${SELF_CONTAINED_PREFIX}"
export NCL_ALREADY_INSTALLED="true"
export PHYCAS_ROOT="${SELF_CONTAINED_BUILD_ROOT}/phycas_trunk"
export PYTHONPATH="${SELF_CONTAINED_PREFIX}:${SELF_CONTAINED_PREFIX}/lib/python2.7/site-packages:${PYTHONPATH}"

