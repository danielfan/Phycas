@echo off
REM
REM Steps to creating a Windows distribution:
REM   1. Build Phycas using bjam
REM   2. Make sure "package_data=windows_package_data" is uncommented in setup.py
REM   3. Change --target-version below to reflect correct Python target version
REM   4. Run this batch file to generate installer (which will be in dist folder)
REM
REM cd ..\..
rmdir /s /q build
rmdir /s /q dist
REM python setup.py sdist --force-manifest
REM copy /Y bin\boost\libs\python\build\boost_python.dll\vc-7_1\release\threading-multi\boost_python.dll phycas\Conversions
python setup.py bdist_wininst --target-version 2.7 --bitmap phycaslogo.bmp --install-script win_shortcuts.py
REM del phycas\__init__.py 
REM svn update
REM cd installer\win
pause
