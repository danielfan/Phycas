#!/bin/sh
set -x
tar xfjv PhycasGUI.app.tar.bz || exit 1
rm -rf ./Phycas.app
mv PhycasGUI.app Phycas.app
cp PhycasGUIexecFromGUI.sh Phycas.app/Contents/MacOS
chmod +x Phycas.app/Contents/MacOS/PhycasGUIexecFromGUI.sh
if test -d phycas_lib
then
	cp -r phycas_lib/* Phycas.app/Contents/Resources/
fi
if test -d phycas_bin
then
	cp -r phycas_bin/* Phycas.app/Contents/MacOS/
fi
if test -d Frameworks
then
	cp -r Frameworks Phycas.app/Contents/Frameworks
fi
cp -r $PHYCAS_ROOT/phycas Phycas.app/Contents/Resources/ || exit 2
find  Phycas.app/Contents/Resources/phycas -name .svn -exec rm -rf {} \; 
sh makeDMG.sh || exit 4

