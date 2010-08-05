#!/bin/sh
sv="$MANUAL_PHYCAS_REVISION_TAG"
if test -z "$sv"
then
	cd "$PHYCAS_ROOT"
	sv=`svnversion -n | sed '/:/d' | sed '/M/d'`
	cd "$PHYCAS_ROOT/installer/mac/iterm"
fi
if test -z "$sv"
then
	echo "Cannot create distribution with out-of-date source"
	echo "MANUAL_PHYCAS_REVISION_TAG is unset and svnversion returned" `svnversion -n ../../../..`
	exit 1
else
	echo "building distro for revision $sv"
fi
set -x

tar xfjv PhycasGUI.app.tar.bz || exit 1
rm -rf ./Phycas.app
mv PhycasGUI.app Phycas.app
cp PhycasGUIexecFromGUI.sh Phycas.app/Contents/MacOS/iTerm4PhycasexecFromGUI.sh
chmod +x Phycas.app/Contents/MacOS/iTerm4PhycasexecFromGUI.sh
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
if test -z "$TARGETTED_PYTHON_VERSION"
then
	echo "TARGETTED_PYTHON_VERSION must be defined to major.minor"
	exit 1
fi 
cp -r "$PHYCAS_ROOT/phycas" Phycas.app/Contents/Resources/ || exit 2
cp -r ./active_phycas_env.sh Phycas.app/Contents/Resources/ || exit 2
cp -r ./run_phycas.sh Phycas.app/Contents/Resources/ || exit 2

cat run_phycas.sh | sed -e "s/PYTHONINTERPRETER=python/PYTHONINTERPRETER=python${TARGETTED_PYTHON_VERSION}/" > Phycas.app/Contents/Resources/run_phycas.sh || exit 2

cp  ./startup.py Phycas.app/Contents/Resources/ || exit 2

if test -z $PY_PACKAGES_TO_BUNDLE
then
	echo not bundling python extras because PY_PACKAGES_TO_BUNDLE is not set!
else
	cp -r $PY_PACKAGES_TO_BUNDLE/setuptools-0.6c11-py2.5.egg Phycas.app/Contents/Resources/ || exit 1
	cp -r $PY_PACKAGES_TO_BUNDLE/IPython Phycas.app/Contents/Resources/ || exit 1
fi


find  Phycas.app/Contents/Resources/phycas -name .svn -exec rm -rf {} \;
cat Phycas.app/Contents/Resources/phycas/__init__.py  | sed "s/PHYCAS_SVN_REVISION_NUMBER_HERE/$sv/g" > t
mv t Phycas.app/Contents/Resources/phycas/__init__.py
sh makeDMG.sh || exit 4

