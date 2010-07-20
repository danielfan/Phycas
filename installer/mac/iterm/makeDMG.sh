#!/bin/sh

set -x
if test -d "$PHYCAS_ROOT"
then
	NAME=Phycas
	if test -f $PHYCAS_ROOT/phycasver.txt
	then
		VERSION=`cat $PHYCAS_ROOT/phycasver.txt`
	else
		echo "Could not find phycasver.txt file and thus could not determine version"
		exit
	fi
	MASTER_DMG=$NAME-$VERSION.dmg	
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
	ditto -rsrc "$PHYCAS_ROOT/LICENSE" "phycasimg/LICENSE" || exit
	ditto -rsrc "README" "phycasimg/README" || exit
	ditto -rsrc "$PHYCAS_ROOT/documentation/users/manual.pdf" "phycasimg/manual.pdf" || exit
	ditto -rsrc "$PHYCAS_ROOT/documentation/users/tutorial/basic.py" "phycasimg/tutorial/basic.py" || exit
	ditto -rsrc "$PHYCAS_ROOT/documentation/users/tutorial/summarize.py" "phycasimg/tutorial/summarize.py" || exit
	ditto -rsrc "$PHYCAS_ROOT/documentation/users/tutorial/partition.py" "phycasimg/tutorial/partition.py" || exit
	ditto -rsrc "$PHYCAS_ROOT/documentation/users/tutorial/steppingstone.py" "phycasimg/tutorial/steppingstone.py" || exit
	ditto -rsrc "$PHYCAS_ROOT/documentation/users/tutorial/polytomy.py" "phycasimg/tutorial/polytomy.py" || exit
	# mkdir "phycasimg/PhycasTutorial"
	# mkdir "phycasimg/PhycasTutorial/scripts"
	# for f in $PHYCAS_ROOT/phycas/Examples/Tutorial/PhycasWHExamples/*nex
	# do 
	# 	ditto -rsrc "$f" "phycasimg/PhycasTutorial" || exit
	# done
	# for f in $PHYCAS_ROOT/phycas/Examples/Tutorial/PhycasWHExamples/scripts/*py
	# do 
	# 	ditto -rsrc "$f" "phycasimg/PhycasTutorial/scripts/" || exit
	# 	ditto -rsrc "$f" "phycasimg/PhycasTutorial/" || exit
	# done
	# ditto -rsrc "$PHYCAS_ROOT/phycas/Examples/Tutorial/phycas_woods_hole_08.pdf" "phycasimg/PhycasTutorial/PhycasTutorial.pdf" || exit
	
	hdiutil create -srcfolder phycasimg $MASTER_DMG
else
	echo "PHYCAS_ROOT is not defined as a directory"
	exit 1
fi
exit 0

