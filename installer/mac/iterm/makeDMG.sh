#!/bin/sh
NAME=Phycas
VERSION=1.0
MASTER_DMG=$NAME-$VERSION.dmg
set -x
if test -d "$PHYCAS_ROOT"
then
	SOURCE_DIR="$PHYCAS_ROOT/installer/mac/iterm"
	cd $SOURCE_DIR || exit
	if test -f $MASTER_DMG
	then
		rm  $MASTER_DMG || exit
	fi
	if test -d ./phycasimg
	then
		rm -r ./phycasimg || exit
	fi
	mkdir phycasimg || exit
	for i in "Phycas.app" 
	do
		ditto -rsrc "./$i" "phycasimg/$i" || exit
	done
	for i in "LICENSE" "readme.txt"
	do
		ditto -rsrc "$PHYCAS_ROOT/$i" "phycasimg/$i" || exit
	done
	hdiutil create -srcfolder phycasimg $MASTER_DMG
else
	echo "PHYCAS_ROOT is not defined as a directory"
	exit 1
fi
exit 0

