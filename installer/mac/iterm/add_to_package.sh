#!/bin/sh
set -x
for i in Development Deployment ; 
do
	cp PhycasGUIexecFromGUI.sh ./build/$i/PhycasGUI.app/Contents/MacOS/ 
	chmod +x ./build/$i/PhycasGUI.app/Contents/MacOS/PhycasGUIexecFromGUI.sh 
	cp -r $PHYCAS_ROOT/phycas  ./build/$i/PhycasGUI.app/Contents/Resources/
done
