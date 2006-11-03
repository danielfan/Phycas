@echo off

REM set PATH=%PATH%;%HOMEDRIVE%\boost_1_33_1\tools\build\jam_src\bin.ntx86
set PHYCAS_ROOT=%HOMEDRIVE%\Synchronized\Projects\ggmods_branch
REM set PYTHON_ROOT=%HOMEDRIVE%\Python24
REM set PYTHON_VERSION=2.4
REM set BOOST_ROOT=$HOMEDRIVE\boost_1_33_1

REM Just in case someone installs VC in a strange place, this will set the necessary
REM environmental variables correctly so that the compiler and its associated DLLs
REM can be located
call "%VS71COMNTOOLS%\..\..\VC7\bin\vcvars32.bat"

if defined PYTHON_ROOT (
	echo PYTHON_ROOT=%PYTHON_ROOT%
) else (
	echo Please define PYTHON_ROOT environmental variable as path to Python24 directory
	echo e.g. set PYTHON_ROOT=C:\Python24
	pause
	exit
)

if defined PYTHON_VERSION (
	echo PYTHON_VERSION=%PYTHON_VERSION%
) else (
	echo Please define PYTHON_VERSION environmental variable as Python major version
	echo e.g. set PYTHON_VERSION=2.4
	pause
	exit
)

if defined BOOST_ROOT (
	echo BOOST_ROOT=%BOOST_ROOT%
) else (
	echo Please define BOOST_ROOT environmental variable as path to boost directory
	echo e.g. set BOOST_ROOT=C:\boost_1_33_0
	pause
	exit
)

REM if defined CIPRES_LIB_ROOT (
REM 	echo CIPRES_LIB_ROOT=%CIPRES_LIB_ROOT%
REM ) else (
REM 	echo Please define CIPRES_LIB_ROOT environmental variable as path to cipres framework directory
REM 	echo e.g. set CIPRES_LIB_ROOT=C:\cipres\framework
REM 	pause
REM 	exit
REM )

if defined PHYCAS_ROOT (
	echo PHYCAS_ROOT=%PHYCAS_ROOT%
) else (
	echo Please define PHYCAS_ROOT environmental variable as path to parent of phycas directory
	echo e.g. set PHYCAS_ROOT=C:\phycasdev
	pause
	exit
)

REM if defined NCL_ROOT (
REM 	echo NCL_ROOT=%NCL_ROOT%
REM ) else (
REM 	echo Please define NCL_ROOT environmental variable as path to parent of ncl directory
REM 	echo e.g. set NCL_ROOT=C:\phycasdev
REM 	pause
REM 	exit
REM )

REM Choose release build and Visual C++ version 7.1 (a.k.a. Visual Studio .NET 2003)
set BUILD=release
echo BUILD=%BUILD%

set TOOLS=vc-7_1
echo TOOLS=%TOOLS%

REM Run bjam, stopping on first error
bjam -q
pause
