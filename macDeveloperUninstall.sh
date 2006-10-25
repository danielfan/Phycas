#!/bin/sh
set -x
sp="$PYTHON_ROOT"/`python  -c 'import distutils.sysconfig;  print distutils.sysconfig.get_python_lib(0,0,prefix="");'`
oldPhypy="$sp/phypy"
lbd=/usr/local/lib/libboost_python.dylib
if test -e "$lbd" ;
then
	sudo rm  "$lbd" || exit
fi

if test -e "$oldPhypy" ;
then
	rm -r "$oldPhypy"
fi

