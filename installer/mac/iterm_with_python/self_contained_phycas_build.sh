#!/bin/sh
export SELF_CONTAINED_BUILD_ROOT=/Users/mholder/Documents/projects/phycas_dev/self_contained
export SELF_CONTAINED_PREFIX="${SELF_CONTAINED_BUILD_ROOT}/PhycasGUI.app/Contents/Resources"
export DYLD_LIBRARY_PATH="$SELF_CONTAINED_PREFIX/lib:$DYLD_LIBRARY_PATH"
export BOOST_ROOT="${SELF_CONTAINED_BUILD_ROOT}/boost_1_43_0"
export BOOST_BUILD_PATH="${BOOST_ROOT}/tools/build/v2"
export PHYCAS_ROOT="${SELF_CONTAINED_BUILD_ROOT}/phycas_trunk"
export OSTYPE="darwin"
export PATH="${BOOST_ROOT}/tools/jam/src/bin.macosxx86_64:$SELF_CONTAINED_PREFIX/bin:$PATH"
export NCL_INSTALL_DIR="${SELF_CONTAINED_PREFIX}"
export NCL_ALREADY_INSTALLED="true"
export PHYCAS_ROOT="${SELF_CONTAINED_BUILD_ROOT}/phycas_trunk"
export PYTHONPATH="${SELF_CONTAINED_PREFIX}:${SELF_CONTAINED_PREFIX}/lib/python2.7/site-packages:${PYTHONPATH}"
