@echo off
REM Assumes cwd is phycasdev\executor and %PHYCAS_ROOT% env var set
REM and cvs.exe is in path
REM
REM **************************************************
echo Get code from CVS...
REM **************************************************
REM cd ..
REM cmd.exe /c cvs update -d
REM cd executor
echo
REM **************************************************
echo Running scomp...
REM **************************************************
cmd.exe /c scomp -out ..\gui\commandLanguage.jar ..\command_archive\xml\phycas_command_language.xsd
cmd.exe /c scomp -out ..\gui\commandState.jar ..\command_archive\xml\command_state.xsd
echo
REM **************************************************
echo Running transform_cmd_xml.py...
REM **************************************************
cd ..
python transform_cmd_xml.py
echo
REM **************************************************
echo Building SocketDebug build in VC...
REM **************************************************
cmd.exe /c devenv /build SocketDebug /project "%PHYCAS_ROOT%\build\VC\phycas\phycas.vcproj" "%PHYCAS_ROOT%\build\VC\VS.sln"
echo
REM **************************************************
echo Copying PhycasSocket.exe to executor directory...
REM **************************************************
cd build
if exist VC\phycas\SocketDebug\phycas.exe (
if exist obsoletePhycasSocket.exe (del obsoletePhycasSocket.exe)
if exist PhycasSocket.exe (ren PhycasSocket.exe obsoletePhycasSocket.exe)
copy VC\phycas\SocketDebug\phycas.exe ..\executor\PhycasSocket.exe
) else (
echo phycas.exe does not exist so it could not be copied
)
echo
REM **************************************************
echo Creating help search database ...
REM **************************************************
cd ..\gui\help
cmd.exe /c jhindexer master
echo
REM **************************************************
echo Running ant...
REM **************************************************
cd ..\..\executor
cmd.exe /c ant -f ..\build\buildGui.xml
cmd.exe /c ant -f ..\build\buildExec.xml
cd ..
echo
REM **************************************************
echo Finished.
REM **************************************************
