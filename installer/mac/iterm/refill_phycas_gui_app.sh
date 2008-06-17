#!/bin/sh
set -x
tar xfjv PhycasGUI.app.tar.bz || exit 1
mv PhycasGUI.app Phycas.app
cp PhycasGUIexecFromGUI.sh Phycas.app/Contents/MacOS
chmod +x Phycas.app/Contents/MacOS/PhycasGUIexecFromGUI.sh
cp -r $PHYCAS_ROOT/phycas Phycas.app/Contents/Resources || exit 2
find  Phycas.app/Contents/Resources/ -name .svn -exec rm -rf {} \; 
sh makeDMG.sh || exit 4
