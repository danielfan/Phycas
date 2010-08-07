@echo off

set PYTHONPATH=%~dp0
set PATH=%~dp0\phycas\Conversions;%PATH%

if "%1"=="" goto noargs
if "%~x1"==".py" goto pydrop

echo The file dropped (%~n1%~x1) does not have the expected .py extension
echo I will try to run it anyway, but if this file does not contain Python
echo code, you will get error messages.
pause

:pydrop 
cd %~p1
python %~n1%~x1
pause
goto :eof

:noargs
python -ic "from phycas import *"
goto :eof
