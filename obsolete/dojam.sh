#!/bin/sh

export PATH=$PATH:$HOME/boost_1_33_1/tools/build/jam_src/bin.linuxx86
export PHYCAS_ROOT=$HOME/Phycas-1.0
export PYTHON_ROOT=/usr
export PYTHON_VERSION=2.4
export BOOST_ROOT=$HOME/boost_1_33_1

if [ -n "$PYTHON_ROOT" ]
then
	echo PYTHON_ROOT=$PYTHON_ROOT
else
	echo Please define PYTHON_ROOT environmental variable as path to Python24 directory
	echo For example:
	echo   if using bash: export PYTHON_ROOT=/Library/Frameworks/Python.framework/Versions/2.4
	echo   if using tcsh: setenv PYTHON_ROOT /Library/Frameworks/Python.framework/Versions/2.4
	exit 0
fi

if [ -n "$PYTHON_VERSION" ]
then
	echo PYTHON_VERSION=$PYTHON_VERSION
else
	echo Please define PYTHON_VERSION environmental variable as Python major version
	echo For example:
	echo   if using bash: export PYTHON_VERSION=2.4
	echo   if using tcsh: setenv PYTHON_VERSION 2.4
	exit 0
fi

if [ -n "$BOOST_ROOT" ]
then
	echo BOOST_ROOT=$BOOST_ROOT
else
	echo Please define BOOST_ROOT environmental variable as path to boost directory
	echo For example:
	echo   if using bash: export BOOST_ROOT=$HOME/boost_1_33_1
	echo   if using tcsh: setenv BOOST_ROOT $HOME/boost_1_33_1
	exit 0
fi

# if [ -n "$CIPRES_LIB_ROOT" ]
# then
# 	echo CIPRES_LIB_ROOT=$CIPRES_LIB_ROOT
# else
# 	echo Please define CIPRES_LIB_ROOT environmental variable as path to cipres framework directory
# 	echo For example:
# 	echo   if using bash: export CIPRES_LIB_ROOT=$HOME/cipres/framework
# 	echo   if using tcsh: setenv CIPRES_LIB_ROOT $HOME/cipres/framework
# 	exit 0
# fi

if [ -n "$PHYCAS_ROOT" ]
then
	echo PHYCAS_ROOT=$PHYCAS_ROOT
else
	echo Please define PHYCAS_ROOT environmental variable as path to parent of phycas directory
	echo For example:
	echo   if using bash: export PHYCAS_ROOT=$HOME/phycasdev
	echo   if using tcsh: setenv PHYCAS_ROOT $HOME/phycasdev
	exit 0
fi

#if [ -n "$NCL_ROOT" ]          
#then
#	echo NCL_ROOT=$NCL_ROOT
#else
#	echo Please define NCL_ROOT environmental variable as path to parent of ncl directory
#	echo For example:
#	echo   if using bash: export NCL_ROOT=$HOME/phycasdev/phypy/src
#	echo   if using tcsh: setenv NCL_ROOT $HOME/phycasdev/phypy/src
#	exit 0
#fi

# export PATH=.:$HOME/local/bin:/usr/local/condor/bin:$PATH
# export DYLD_LIBRARY_PATH=$HOME/phycasdev/phypy/bin/boost/libs/python/build/libboost_python.dylib/darwin/release/shared-linkable-true
# export CVS_RSH=/usr/bin/ssh
# export PYTHONPATH=$HOME/site-packages/numarray:$HOME/phycasdev/phypy/phypy:$HOME/phycasdev/phypy/apps:$HOME/phycasdev/python

export BUILD=release

if [ "$OSTYPE" = "darwin" ]
then          
	export TOOLS=darwin
else
	export TOOLS=gcc
fi

bjam -q
