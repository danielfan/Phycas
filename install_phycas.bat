@echo off
echo Setting environmental variable PYTHONPATH to...
echo %CD%
setx PYTHONPATH "%CD%" > NUL
echo.
echo Adding this to beginning of environmental variable PATH...
echo %CD%\phycas\Conversions
setx PATH "%CD%\phycas\Conversions;%PATH%" > NUL
echo.
echo Phycas can now be invoked from any command window.
pause
