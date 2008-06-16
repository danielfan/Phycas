#!/bin/sh
mac_os_dir=`dirname $0`
contents_os_dir=`dirname $mac_os_dir`
resources_dir="$contents_os_dir/Resources"
phycas_dir="$resources_dir/phycas"
if test -d $phycas_dir
then
	if test -z $PYTHONPATH
	then
		PYTHONPATH="$resources_dir"
	else
		PYTHONPATH="$PYTHONPATH:$resources_dir"
	fi
	if test -z $DYLD_LIBRARY_PATH
	then
		DYLD_LIBRARY_PATH="$phycas_dir/Conversions"
	else
		DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$phycas_dir/Conversions"
	fi
else
	echo "$0: $phycas_dir does not exist!"
	exit 1
fi
exit 0
