#!/bin sh
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


if ! test -f Python-2.6.5.tar.bz2
then
    wget http://www.python.org/ftp/python/2.6.5/Python-2.6.5.tar.bz2 || exit
fi


if ! test  -f boost_1_43_0.tar.bz2
then
    wget http://sourceforge.net/projects/boost/files/boost/1.43.0/boost_1_43_0.tar.bz2/download || exit
fi

mkdir -p "${SELF_CONTAINED_BUILD_ROOT}"
cd "${SELF_CONTAINED_BUILD_ROOT}" || exit 1


if ! test -d trunk_phycas
then
    svn co https://phycas.svn.sourceforge.net/svnroot/phycas/trunk trunk_phycas || exit
fi
if ! test -d PhycasGUI.app
then
    tar xfvj trunk_phycas/installer/mac/iterm/PhycasGUI.app.tar.bz  || exit
fi


if ! test -d readline-6.1
then
    cp "${script_dir}/readline-6.1.tar" .
    tar xfv readline-6.1.tar || exit 1
fi
cd readline-6.1  || exit 1
CFLAGS="$CFLAGS" LDFLAGS="$LDFLAGS" ./configure --prefix=$SELF_CONTAINED_PREFIX || exit 1
make  -j4 || exit 1
make check  || exit 1
make install  || exit 1
cd ..  || exit 1


if ! test -d Python-2.6.5
then
    tar xfvj "${script_dir}/Python-2.6.5.tar.bz2"  || exit
fi
cd  Python-2.6.5  || exit
export EXTRA_CFLAGS="$CFLAGS -isysroot /Developer/SDKs/MacOSX10.5.sdk -I${SELF_CONTAINED_PREFIX}/include -I${SELF_CONTAINED_PREFIX}/include/readline " LDFLAGS="$LDFLAGS -L${SELF_CONTAINED_PREFIX}/lib" 
./configure --prefix="${SELF_CONTAINED_PREFIX}" --enable-shared

make -j4 || exit
make install  || exit
cd .. || exit


if ! test -d nclv2.1
then
    svn co https://ncl.svn.sourceforge.net/svnroot/ncl/branches/v2.1 nclv2.1  || exit 1
fi
cd nclv2.1 || exit 1
sh bootstrap.sh || exit 1
LDFLAGS="$LDFLAGS" CXXFLAGS="$CXXFLAGS -isysroot /Developer/SDKs/MacOSX10.5.sdk  -Wreturn-type -Woverloaded-virtual -Wformat -Wmissing-braces -Wparentheses -Wswitch -Wunused-function -Wunused-label -Wunused-parameter -Wunused-variable -Wunused-value -Wunknown-pragmas -pedantic -Wshadow -Wfour-char-constants -Wsign-compare -Wnewline-eof -Wall -Wreturn-type -Wunused -Wredundant-decls -Wcast-align -Wcomment -Wextra" ./configure --prefix=${SELF_CONTAINED_PREFIX} --disable-static || exit 1
make -j4 || exit 1
make check || exit 1
make install || exit 1
make installcheck || exit 1
rm ${SELF_CONTAINED_PREFIX}/bin/NCLconverter 
rm ${SELF_CONTAINED_PREFIX}/bin/NEXUSnormalizer 
rm ${SELF_CONTAINED_PREFIX}/bin/NEXUSvalidator
cd ..


