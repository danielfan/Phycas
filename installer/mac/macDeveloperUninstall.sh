#!/bin/sh
set -x
sp="$PYTHON_ROOT"/`python  -c 'import distutils.sysconfig;  print distutils.sysconfig.get_python_lib(0,0,prefix="");'`
oldPhypy="$sp/phypy"
oldPhypyEgg="$sp/Phycas-1.0-py${PYTHON_VERSION}.egg-info"
lbd=/usr/local/lib/libboost_python.dylib
if test -e "$lbd" ;
then
	sudo rm  "$lbd" || exit
fi

if test -e "$oldPhypy" ;
then
	rm -rf "$oldPhypy" || exit
fi

if test -e "$oldPhypyEgg" ;
then
	rm -rf "$oldPhypyEgg" || exit
fi

