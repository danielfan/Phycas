#!/bin/sh
set -x
tar xfjv PhycasGUI.app.tar.bz || exit 1
cp -r $PHYCAS_ROOT/phycas PhycasGUI.app/Contents/Resources || exit 2
find  PhycasGUI.app/Contents/Resources/ -name .svn -exec rm -rf {} \; 
sh makeDMG.sh || exit 4
