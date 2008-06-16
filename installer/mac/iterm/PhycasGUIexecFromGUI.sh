#!/bin/sh
env
mac_os_dir=`dirname $0`
contents_os_dir=`dirname $mac_os_dir`
phycas_dir=`$contents_os_dir/Resources/phycas"
if test -d $phycas_dir
then
else
	echo "$0: $phycas_dir does not exist!"
	exit 1
fi
exit 0
