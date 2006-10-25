#!/bin/sh
# Steps to creating a  distribution:
#   1. invoke  phypy/dojam.py
#   2. Copy libboost_python.dylib to phypy/phypy/Conversions 
#   3. Make sure "package_data=windows_package_data" is uncommented in setup.py
#   4. Change --target-version below to reflect correct Python target version
#   5. Run this batch file to generate installer (which will be in dist folder)
set -x
rm -rf build
cd phypy
python dojam.py || exit
cd ..
cp phypy/bin/boost/libs/python/build/libboost_python.dylib/darwin/release/shared-linkable-true/libboost_python.dylib phypy/phypy/Conversions
if test -e "$PYTHON_ROOT/bin/bdist_mpkg" ;
then
    python setup.py build --compiler=fake
    "$PYTHON_ROOT/bin/bdist_mpkg" --skip-build -k
else
    echo "bdist_mpkp is required.  see http://cheeseshop.python.org/pypi/bdist_mpkg/"
    exit 1
fi
