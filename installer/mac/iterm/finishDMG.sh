#!/bin/sh
NAME=Phycas
VERSION=1.0
MASTER_DMG=$NAME-$VERSION.dmg
set -x
if test -d "$PHYCAS_ROOT"
then
	cd $PHYCAS_ROOT/installer/mac/diskimage/
	rm -f "$MASTER_DMG"
	hdiutil convert "$WC_DMG" -quiet -format UDZO -imagekey zlib-level=9 -o "$@"
	#rm -rf $WC_DIR
	mv "$WC_DMG" "$MASTER_DMG"
	mv "$MASTER_DMG" "$SOURCE_DIR/"
else
	echo "PHYCAS_ROOT is not defined as a directory"
	exit 1
fi
exit 0
