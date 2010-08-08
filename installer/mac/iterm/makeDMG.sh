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
	mkdir "phycasimg/tutorial"
	mkdir "phycasimg/tutorial/ex1basic"
	mkdir "phycasimg/tutorial/ex2steppingstone"
	mkdir "phycasimg/tutorial/ex3-4polytomyPriors"
	for f in ex1basic/basic.py ex1basic/green.nex ex2steppingstone/murphy29.rag1rag2.nex ex2steppingstone/setUpByCodonPosPartition.py ex2steppingstone/setUpByGeneAndCodonPosPartition.py ex2steppingstone/steppingstone1.py ex2steppingstone/steppingstone2.py ex2steppingstone/steppingstone3.py ex2steppingstone/steppingstoneDoMCMC.py ex2steppingstone/steppingstoneReadData.py ex3-4polytomyPriors/IntExtPrior.py ex3-4polytomyPriors/NoPolytomy.py ex3-4polytomyPriors/Polytomy.py ex3-4polytomyPriors/ShoupLewis.nex
	do
		ditto -rsrc "$PHYCAS_ROOT/documentation/users/tutorial/$i" "phycasimg/tutorial/$i" || exit
	done

	hdiutil create -srcfolder phycasimg $MASTER_DMG
else
	echo "PHYCAS_ROOT is not defined as a directory"
	exit 1
fi
exit 0

