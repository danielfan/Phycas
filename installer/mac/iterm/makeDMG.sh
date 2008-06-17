#!/bin/sh
NAME=Phycas
VERSION=1.0
MASTER_DMG=$NAME-$VERSION.dmg
set -x
files="PhycasGUI.app" 
if test -d "$PHYCAS_ROOT"
then
	SOURCE_DIR="$PHYCAS_ROOT/installer/mac/iterm"
	cd $PHYCAS_ROOT/installer/mac/diskimage/
	TEMPLATE_DMG="template.dmg"
	bunzip2 -k ${TEMPLATE_DMG}.bz2
	WC_DMG=wapp.dmg
	WC_DIR=wapp
	cp $TEMPLATE_DMG $WC_DMG
	if test ./$WC_DIR
	then
		rm -rf ./$WC_DIR
	fi
	mkdir -p ${WC_DIR}
	hdiutil attach "$WC_DMG" -noautoopen -quiet -mountpoint "$WC_DIR"
	for i in $files
	do
		rm -rf "$WC_DIR/$i"; \
		ditto -rsrc "$SOURCE_DIR/$i" "$WC_DIR/$i"; \
	done
	WC_DEV=`hdiutil info | grep "$WC_DIR" | grep "Apple_HFS" | awk '{print $1}'` 
	if test -z $WC_DEV
	then
		WC_DEV=`hdiutil info | grep "$WC_DIR" | grep '^/' | awk '{print $1}'` 
	fi
	if test -z $WC_DEV
	then
		echo "Could not unmount. Run 'hdiutil info' "
		exit 2
	else
		hdiutil detach $WC_DEV -quiet -force
		rm -f "$MASTER_DMG"
		hdiutil convert "$WC_DMG" -quiet -format UDZO -imagekey zlib-level=9 -o "$@"
		#rm -rf $WC_DIR
		mv "$WC_DMG" "$MASTER_DMG"
		mv "$MASTER_DMG" "$SOURCE_DIR/"
	fi
else
	echo "PHYCAS_ROOT is not defined as a directory"
	exit 1
fi
exit 0
