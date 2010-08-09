copy /b MANIFEST.common+MANIFEST.windows MANIFEST.in
python setup.py sdist --formats=zip
pause
