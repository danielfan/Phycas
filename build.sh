#!/bin/sh

export BOOST_ROOT="$HOME/boost_1_36_0"
export PATH="${PATH}:$BOOST_ROOT/tools/jam/src/bin.linuxx86"
export BOOST_BUILD_PATH="${BOOST_ROOT}/tools/build/v2"
export NCL_INSTALL_DIR="$HOME/ncldev/branches/v2.1"
export PHYCAS_ROOT="$HOME/phycasdev/trunk"
bjam release
