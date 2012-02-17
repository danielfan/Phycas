#!/bin/sh

export OSTYPE="darwin"
export BOOST_ROOT="/Users/Daniel/Desktop/Phycas/boost_1_48_0"
export PATH="${PATH}:$BOOST_ROOT/tools/build/v2/engine/bin.macosxx86_64"
export BOOST_BUILD_PATH="${BOOST_ROOT}/tools/build/v2"
export PHYCAS_ROOT="/Users/Daniel/Desktop/Phycas/Phycas"
export NCL_ALREADY_INSTALLED=1
export NCL_INSTALL_DIR="/Users/Daniel/Desktop/Phycas"
bjam -q release
