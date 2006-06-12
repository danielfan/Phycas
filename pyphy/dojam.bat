@echo off

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

if defined CIPRES_LIB_ROOT (
	echo CIPRES_LIB_ROOT=%CIPRES_LIB_ROOT%
) else (
	echo Please define CIPRES_LIB_ROOT environmental variable as path to cipres framework directory
	echo e.g. set CIPRES_LIB_ROOT=C:\cipres\framework
	pause
	exit
)

if defined PHYCAS_ROOT (
	echo PHYCAS_ROOT=%PHYCAS_ROOT%
) else (
	echo Please define PHYCAS_ROOT environmental variable as path to parent of phycas directory
	echo e.g. set PHYCAS_ROOT=C:\phycasdev
	pause
	exit
)

if defined NCL_ROOT (
	echo NCL_ROOT=%NCL_ROOT%
) else (
	echo Please define NCL_ROOT environmental variable as path to parent of ncl directory
	echo e.g. set NCL_ROOT=C:\phycasdev
	pause
	exit
)

REM Choose release build and Visual C++ version 7.1 (a.k.a. Visual Studio .NET 2003)
set BUILD=release
echo BUILD=%BUILD%

set TOOLS=vc-7_1
echo TOOLS=%TOOLS%

REM Run bjam, stopping on first error
bjam -q
pause
