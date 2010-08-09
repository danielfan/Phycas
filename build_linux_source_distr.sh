#!/bin/sh 

rm -rf dist
rm -f MANIFEST.in
cat MANIFEST.common > MANIFEST.in
cat MANIFEST.linux >> MANIFEST.in
python setup.py sdist --force-manifest --formats=gztar
