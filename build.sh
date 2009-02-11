#!/bin/sh

export OSTYPE="darwin"
export BOOST_ROOT="$HOME/boost_1_37_0"
export PATH="${PATH}:$BOOST_ROOT/tools/jam/src/bin.macosxx86"
export BOOST_BUILD_PATH="${BOOST_ROOT}/tools/build/v2"
export NCL_INSTALL_DIR="../ncl-2.1.382"
export PHYCAS_ROOT="."
bjam release
