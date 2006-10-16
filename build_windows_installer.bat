@echo off
REM
REM Steps to creating a Windows distribution:
REM   1. Build both boost_python.dll and Phycas using dojam.bat
REM   2. Copy boost_python.dll to phycasdev/phypy/phypy/Conversions 
REM   3. Make sure "package_data=windows_package_data" is uncommented in setup.py
REM   4. Change --target-version below to reflect correct Python target version
REM   5. Run this batch file to generate installer (which will be in dist folder)
REM
rmdir /s /q build
rmdir /s /q dist
REM python setup.py sdist --force-manifest
copy /Y phypy\bin\boost\libs\python\build\boost_python.dll\vc-7_1\release\threading-multi\boost_python.dll phypy\phypy\Conversions
python setup.py bdist_wininst --target-version 2.4 --bitmap phycaslogo.bmp --install-script win_shortcuts.py
pause
