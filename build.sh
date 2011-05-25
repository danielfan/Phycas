#!/bin/sh

export OSTYPE="darwin"
export BOOST_ROOT="$HOME/tools/boost_1_41_0"
export PATH="${PATH}:$BOOST_ROOT/tools/jam/src/bin.macosxx86_64"
export BOOST_BUILD_PATH="${BOOST_ROOT}/tools/build/v2"
export NCL_INSTALL_DIR="../ncl-svn/branches/v2.2"
export PHYCAS_ROOT="."
bjam release
