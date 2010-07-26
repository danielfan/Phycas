#!/bin sh
set -x 
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
    wget ftp://ftp.cwru.edu/pub/bash/readline-6.1.tar.gz || exit
    tar xfvz readline-6.1.tar.gz || exit
fi
cd readline-6.1  || exit
./configure --prefix=$SELF_CONTAINED_PREFIX || exit
make && make check && make install  || exit
cd ..  || exit

if ! test -d Python-2.7
then
    wget http://www.python.org/ftp/python/2.7/Python-2.7.tar.bz2 || exit
    tar xfvj Python-2.7.tar.bz2  || exit
fi
cd  Python-2.7  || exit
CCFLAGS="-L ${SELF_CONTAINED_PREFIX}/lib -I ${SELF_CONTAINED_PREFIX}/include ${SELF_CONTAINED_PREFIX}/include/readline" ./configure --prefix="${SELF_CONTAINED_PREFIX}" --enable-framework=/Users/mholder/Documents/projects/phycas_dev/self_contained/PhycasGUI.app/Contents/Framework || exit
make || exit
make install  || exit
cd .. || exit

if ! test  -d boost_1_43_0
then
    wget http://sourceforge.net/projects/boost/files/boost/1.43.0/boost_1_43_0.tar.bz2/download || exit
    tar xfjv boost_1_43_0.tar.bz2  || exit
fi
cd boost_1_43_0/tools/jam/src || exit
./build.sh || exit
cd ../../../.. || exit


if ! test -d nclv2.1
then
    svn co https://ncl.svn.sourceforge.net/svnroot/ncl/branches/v2.1 nclv2.1 || exit
fi
cd nclv2.1
sh bootstrap.sh
CXXFLAGS="-Wreturn-type -Woverloaded-virtual -Wformat -Wmissing-braces -Wparentheses -Wswitch -Wunused-function -Wunused-label -Wunused-parameter -Wunused-variable -Wunused-value -Wunknown-pragmas -pedantic -Wshadow -Wfour-char-constants -Wsign-compare -Wnewline-eof -Wall -Wreturn-type -Wunused -Wredundant-decls -Wcast-align -Wcomment -Wextra" ./configure --prefix=${SELF_CONTAINED_PREFIX} --disable-static || exit
make || exit
make check || exit
make install || exit
make installcheck || exit
cd .. || exit

cd trunk_phycas || exit
bjam release -q || exit
cp -r phycas  $SELF_CONTAINED_PREFIX/phycas
cd ..


if ! test -f setuptools-0.6c11-py2.7.egg
then
    wget http://pypi.python.org/packages/2.7/s/setuptools/setuptools-0.6c11-py2.7.egg#md5=fe1f997bc722265116870bc7919059ea
fi
sh setuptools-0.6c11-py2.7.egg --prefix=$SELF_CONTAINED_PREFIX



if ! test -d ipython-0.10
then
    wget http://ipython.scipy.org/dist/0.10/ipython-0.10.tar.gz
    tar xfvz ipython-0.10.tar.gz
fi
cd ipython-0.10 || exit
python setup.py install --prefix=$SELF_CONTAINED_PREFIX || exit
cd .. || exit

cp  trunk_phycas/installer/mac/iterm_with_python/iTerm4PhycasexecFromGUI.sh PhycasGUI.app/Contents/MacOS/iTerm4PhycasexecFromGUI.sh || exit
