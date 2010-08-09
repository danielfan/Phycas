#!/bin/bash

# This removes dynamic link libraries from previous builds
./remove_dylibs.sh

# Specify "linux" or "darwin" (for macos)
export OSTYPE="linux"

# The path to your boost installation (download from http://www.boost.org/) 
export BOOST_ROOT="$HOME/boost_1_42_0"

# Modify PATH so that bjam executable can be found
# First, create bjam executable by running $BOOST_ROOT/tools/jam/src/build.sh
# The build.sh script will create a macosxx86_64 folder on mac snow leopard
# and a linuxx86 folder on linux. Be sure final part of path below is correct
# for your platform.
export PATH="${PATH}:$BOOST_ROOT/tools/jam/src/bin.linuxx86"

# This is needed for bjam to find its way
export BOOST_BUILD_PATH="${BOOST_ROOT}/tools/build/v2"

# Provide path to Phycas source code
export PHYCAS_ROOT="$HOME/Phycas-1.2.0"

# Provide path to preinstalled Nexus Class Library 
# (download from http://sourceforge.net/projects/ncl/)
export NCL_ALREADY_INSTALLED=1
export NCL_INSTALL_DIR="/usr/local"

# This command initiates the build
bjam release
