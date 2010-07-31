#!/bin sh
set -x
if ! test -f concrete_phycas_build_env.sh
then 
    echo "Running env.sh to create concrete_phycas_build_env.sh"
    sh env.sh || exit 1
fi
if test -f concrete_phycas_build_env.sh
then 
    source concrete_phycas_build_env.sh
else
    echo "concrete_phycas_build_env.sh not found!"
    exit 1
fi

if test -z ${SELF_CONTAINED_BUILD_ROOT}
then
    echo "SELF_CONTAINED_BUILD_ROOT must be defined"
    exit 1
fi
if test -z $SELF_CONTAINED_PREFIX
then
    echo "SELF_CONTAINED_PREFIX must be defined"
    exit 1
fi

rel_script_dir=`dirname $0`
if ! test -d $rel_script_dir
then
    echo "$rel_script_dir not found"
    exit 1
fi
script_dir=$(python -c "import os; print os.path.abspath('$rel_script_dir')")



################################################################################
echo "Downloading dependencies (if needed)."
if ! test -f "readline-6.1.tar.gz"
then 
    wget ftp://ftp.cwru.edu/pub/bash/readline-6.1.tar.gz || exit 1
fi

if ! test -f Python-2.6.5.tar.bz2
then
    wget http://www.python.org/ftp/python/2.6.5/Python-2.6.5.tar.bz2 || exit 1
fi


if ! test  -f boost_1_43_0.tar.bz2
then
    wget http://sourceforge.net/projects/boost/files/boost/1.43.0/boost_1_43_0.tar.bz2/download || exit 1
fi



if ! test -f ipython-0.10.tar.gz
then
    wget http://ipython.scipy.org/dist/0.10/ipython-0.10.tar.gz || exit 1
fi

################################################################################

mkdir -p "${SELF_CONTAINED_BUILD_ROOT}"
cd "${SELF_CONTAINED_BUILD_ROOT}" || exit 1

if ! test -d trunk_phycas
then
    echo "Checking out Phycas"

    svn co https://phycas.svn.sourceforge.net/svnroot/phycas/trunk trunk_phycas || exit
fi
if ! test -d PhycasGUI.app
then
    echo "Unpacking the PhycasGUI.app"
    tar xfvj trunk_phycas/installer/mac/iterm/PhycasGUI.app.tar.bz  || exit
fi


if ! test -f "$SELF_CONTAINED_PREFIX/lib/libreadline.dylib"
then
    echo "Installing readline into the package"
    if ! test -d readline-6.1
    then
        tar xfvz "${script_dir}/readline-6.1.tar.gz" || exit 1
    fi
    cd readline-6.1  || exit 1
    CFLAGS="$CFLAGS" LDFLAGS="$LDFLAGS" ./configure --prefix="$SELF_CONTAINED_PREFIX" || exit 1
    make  -j4 || exit 1
    make check  || exit 1
    make install  || exit 1
    cd ..  || exit 1
fi

if ! test -f "$SELF_CONTAINED_PREFIX/bin/python"
then
    echo "Installing Pytho into the package."
    echo "We have to use EXTRA_CFLAGS and --enable-shared"
    if ! test -d Python-2.6.5
    then
        tar xfvj "${script_dir}/Python-2.6.5.tar.bz2"  || exit
    fi
    cd  Python-2.6.5  || exit
    export EXTRA_CFLAGS="$CFLAGS -isysroot /Developer/SDKs/MacOSX10.5.sdk -I${SELF_CONTAINED_PREFIX}/include -I${SELF_CONTAINED_PREFIX}/include/readline " LDFLAGS="$LDFLAGS -L${SELF_CONTAINED_PREFIX}/lib" 
    ./configure --prefix="${SELF_CONTAINED_PREFIX}" --enable-shared || exit 1

    # the following 3 lines is a hack to force hashlib and md5 lib to be built
    cat setup.py  | sed -E "s/#print 'openssl_ver.*/openssl_ver=0 # MTH-hack/" > t
    mv setup.py setup.py.BAK
    mv t setup.py
    
    
    make -j4 || exit
    make install  || exit
    cd .. || exit
fi

if ! test -f "$SELF_CONTAINED_PREFIX/lib/ncl/libncl.dylib"
then
    echo "installing NCL"
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
fi

if ! test -f "$BOOST_ROOT/tools/jam/src/bin.macosxx86_64/bjam"
then
    if ! test  -d boost_1_43_0
    then
        tar xfvj "${script_dir}/boost_1_43_0.tar.bz2"  || exit
    fi
    cd $BOOST_ROOT || exit
    #cp "${script_dir}/project-config.jam.bak" project-config.jam || exit 1
    sh bootstrap.sh --with-toolset=darwin --with-libraries=python,thread "--with-python=${SELF_CONTAINED_PREFIX}/bin/python" "--with-python-root=${SELF_CONTAINED_PREFIX}" --with-python-version=2.6 || exit 1
    ./bjam release cxxflags="-m32 -arch i386  -isysroot /Developer/SDKs/MacOSX10.5.sdk" linkflags="-m32 -arch i386" -q -d2 --toolset=darwin
    #./bjam release cxxflags="-m32 -arch i386  -isysroot /Developer/SDKs/MacOSX10.5.sdk -I${SELF_CONTAINED_PREFIX}/include/python2.6" linkflags="-m32 -arch i386 -L${SELF_CONTAINED_PREFIX}/lib -lpython2.6" -q -d2 --toolset=darwin || exit 1
    cd ./tools/jam/src || exit 1
    ./build.sh || exit 1
    cd ../../../.. || exit 1
fi


if ! test -d "$SELF_CONTAINED_PREFIX/lib/python2.6/site-packages/phycas"
then
    cd trunk_phycas || exit 1
    bjam release  cxxflags="-m32 -arch i386  -isysroot /Developer/SDKs/MacOSX10.5.sdk -g -include phycas/src/std_force_include.hpp $CXXFLAGS -I$NCL_INSTALL_DIR/include -I. -I${SELF_CONTAINED_PREFIX}/include/python2.6 " cflags="-m32 -arch i386  -isysroot /Developer/SDKs/MacOSX10.5.sdk -g -include phycas/src/std_force_include.hpp $CXXFLAGS -I$NCL_INSTALL_DIR/include -I. -I${SELF_CONTAINED_PREFIX}/include/python2.6 " linkflags="-m32 -arch i386 -L$NCL_INSTALL_DIR/lib/ncl -lncl -L${SELF_CONTAINED_PREFIX}/lib -lpython2.6" -q -d3  || exit 1
    
    for i in Conversions/_ConversionsExt.so DataMatrix/_DataMatrixExt.so Likelihood/_LikelihoodExt.so Phylogeny/_PhylogenyExt.so ProbDist/_ProbDistExt.so ReadNexus/_ReadNexusExt.so
    do
        if ! test -f phycas/$i
        then
            echo "$i was not built"
            exit 1
        fi
    done
    
    cp -r phycas  "$SELF_CONTAINED_PREFIX/lib/python2.6/site-packages/phycas" || exit 1
    cp "installer/mac/iterm_with_python/iTerm4PhycasexecFromGUI.sh" "${SELF_CONTAINED_CONTENTS_DIR}/MacOS"
    cd ..
fi



if ! test -f "$SELF_CONTAINED_PREFIX/bin/easy_install"
then
    if ! test -f setuptools-0.6c11-py2.6.egg
    then
        wget 'http://pypi.python.org/packages/2.6/s/setuptools/setuptools-0.6c11-py2.6.egg#md5=bfa92100bd772d5a213eedd356d64086'
    fi
    sh setuptools-0.6c11-py2.6.egg --prefix=$SELF_CONTAINED_PREFIX
fi


if ! test -f "$SELF_CONTAINED_PREFIX/bin/ipython"
then
    if ! test -d ipython-0.10
    then
        tar xfvj "${script_dir}/ipython-0.10.tar.gz"
    fi
    cd ipython-0.10 || exit
    python setup.py install --prefix=$SELF_CONTAINED_PREFIX || exit
    cd .. || exit
fi

# CLEANUP
grep -r "#!$SELF_CONTAINED_CONTENTS_DIR/Resources" . | sed -E 's/(\.\/.*):#!\/.*/\1/' > with_full_path_shebang
escaped_c=$(echo $SELF_CONTAINED_CONTENTS_DIR | sed -e 's/\//\\\//g')
for i in `cat with_full_path_shebang`
do
    sed -E "s/#!${escaped_c}\/Resources\/bin\//#!\/usr\/bin\/env /" $i > y
    diff y $i
    cp y $i
done


if test -z "$SELF_CONTAINED_PREFIX"
then 
    echo "SELF_CONTAINED_PREFIX must be defined"
else
    if test -d "$SELF_CONTAINED_PREFIX/lib/python2.6/site-packages/phycas"
    then
        find "$SELF_CONTAINED_PREFIX/lib/python2.6/site-packages/phycas" -name .svn -exec rm -rf {} \;
        rm -r "$SELF_CONTAINED_PREFIX/lib/python2.6/site-packages/phycas/src"
    
    fi
    find "$SELF_CONTAINED_PREFIX" -name "*.pyc" -exec rm {} \;
fi
