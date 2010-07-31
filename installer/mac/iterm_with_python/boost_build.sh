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


