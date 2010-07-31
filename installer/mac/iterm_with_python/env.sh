#!/bin/sh
export PATH=/usr/local/bin:/bin:/sbin:/usr/bin:/usr/sbin
echo "Setting the PATH to a *very* limited. $PATH to avoid us picking up a tool in an odd place"
unset DYLD_LIBRARY_PATH
unset PYTHONPATH

rel_script_dir=`dirname $0`
if ! test -d $rel_script_dir
then
    echo "$rel_script_dir not found"
    exit 1
fi

export SELF_CONTAINED_BUILD_ROOT=$(python -c "import os; print os.path.abspath('$script_dir')")/build
export SELF_CONTAINED_CONTENTS_DIR="${SELF_CONTAINED_BUILD_ROOT}/PhycasGUI.app/Contents"
export SELF_CONTAINED_PREFIX="${SELF_CONTAINED_CONTENTS_DIR}/Resources"


# downstream we assume 2.6, so don't think that changing this will work!!!
pyversion=2.6

export DYLD_LIBRARY_PATH="$SELF_CONTAINED_PREFIX/lib:$SELF_CONTAINED_PREFIX/lib/ncl:$SELF_CONTAINED_PREFIX/lib/python${SELF_CONTAINED_PREFIX}/site-packages/phycas:${DYLD_LIBRARY_PATH}"
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
export PYTHONPATH="${SELF_CONTAINED_PREFIX}:${SELF_CONTAINED_PREFIX}/lib/python${SELF_CONTAINED_PREFIX}/site-packages:${PYTHONPATH}"

echo "Writing concrete_phycas_build_env.sh for future reference"
echo '#!/bin/sh' > concrete_phycas_build_env.sh
for p in SELF_CONTAINED_BUILD_ROOT SELF_CONTAINED_CONTENTS_DIR SELF_CONTAINED_PREFIX DYLD_LIBRARY_PATH CXXFLAGS CFLAGS LDFLAGS MACOSX_DEPLOYMENT_TARGET BOOST_ROOT BOOST_BUILD_PATH PHYCAS_ROOT OSTYPE PATH NCL_INSTALL_DIR NCL_ALREADY_INSTALLED PHYCAS_ROOT PYTHONPATH
do
    v=$(echo \$$p)
    echo "echo export $p=\\\"$v\\\"" > .dummy_shell.sh
    sh .dummy_shell.sh >> concrete_phycas_build_env.sh
    rm .dummy_shell.sh
done
