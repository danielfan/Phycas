#!/bin/sh 

rm -rf dist
python setup.py sdist --force-manifest --formats=gztar
pause
