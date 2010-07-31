#!/bin/sh
set -x 
script_dir=$(python -c "import os; print os.path.split(os.path.abspath(\"$0\"))[0]")
echo $script_dir
if test -z ${SELF_CONTAINED_BUILD_ROOT}
then
    echo "SELF_CONTAINED_BUILD_ROOT must be defined"
    exit 1
fi

if test -z $SELF_CONTAINED_PREFIX
then
    echo SELF_CONTAINED_PREFIX must be defined
    exit 1
fi

mkdir -p "${SELF_CONTAINED_BUILD_ROOT}"
cd "${SELF_CONTAINED_BUILD_ROOT}" || exit 1









cd trunk_phycas || exit
bjam release  cxxflags="-m32 -arch i386  -isysroot /Developer/SDKs/MacOSX10.5.sdk -g -include phycas/src/std_force_include.hpp $CXXFLAGS -I$NCL_INSTALL_DIR/include -I. -I${SELF_CONTAINED_PREFIX}/include/python2.6 " cflags="-m32 -arch i386  -isysroot /Developer/SDKs/MacOSX10.5.sdk -g -include phycas/src/std_force_include.hpp $CXXFLAGS -I$NCL_INSTALL_DIR/include -I. -I${SELF_CONTAINED_PREFIX}/include/python2.6 " linkflags="-m32 -arch i386 -L$NCL_INSTALL_DIR/lib/ncl -lncl -L${SELF_CONTAINED_PREFIX}/lib -lpython2.6" -q -d3  || exit 1

for i in Conversions/_ConversionsExt.so DataMatrix/_DataMatrixExt.so Likelihood/_LikelihoodExt.so Phylogeny/_PhylogenyExt.so ProbDist/_ProbDistExt.so ReadNexus/_ReadNexusExt.so
do
	if ! test -f phycas/$i
	then
		echo "$i was not built"
		exit 1
	fi
done

cp -r phycas  $SELF_CONTAINED_PREFIX/lib/python2.6/site-packages/phycas
cd ..

