@echo off
if exist runall_diffs.txt del runall_diffs.txt

set TESTDIR=%CD%
cd ..\Examples

echo ****************************
echo *** Running ExplorePrior ***
echo ****************************
echo **************************** >> %TESTDIR%\runall_diffs.txt
echo *** Running ExplorePrior *** >> %TESTDIR%\runall_diffs.txt
echo **************************** >> %TESTDIR%\runall_diffs.txt
echo. >> %TESTDIR%\runall_diffs.txt
cd ExplorePrior
start /low /b /wait python ExplorePrior.py
fc nodata.nex.p reference_output\nodata.nex.p >> %TESTDIR%\runall_diffs.txt 
if errorlevel 1 (goto abort)
fc nodata.nex.t reference_output\nodata.nex.t >> %TESTDIR%\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

echo.
echo ***************************
echo *** Running FixedParams ***
echo ***************************
echo. >> runall_diffs.txt
echo *************************** >> %TESTDIR%\runall_diffs.txt
echo *** Running FixedParams *** >> %TESTDIR%\runall_diffs.txt
echo *************************** >> %TESTDIR%\runall_diffs.txt
echo. >> %TESTDIR%\runall_diffs.txt
cd FixedParams
start /low /b /wait python FixedParams.py
fc params.p reference_output\params.p >> %TESTDIR%\runall_diffs.txt
if errorlevel 1 (goto abort)
fc trees.t reference_output\trees.t >> %TESTDIR%\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

echo.
echo ****************************
echo *** Running GelfandGhosh ***
echo ****************************
echo. >> runall_diffs.txt
echo **************************** >> %TESTDIR%\runall_diffs.txt
echo *** Running GelfandGhosh *** >> %TESTDIR%\runall_diffs.txt
echo **************************** >> %TESTDIR%\runall_diffs.txt
echo. >> %TESTDIR%\runall_diffs.txt
cd GelfandGhosh
start /low /b /wait python GelfandGhosh.py
fc ggout.txt reference_output\ggout.txt >> %TESTDIR%\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

echo.
echo ******************************
echo *** Running LikelihoodTest ***
echo ******************************
echo. >> runall_diffs.txt
echo ****************************** >> %TESTDIR%\runall_diffs.txt
echo *** Running LikelihoodTest *** >> %TESTDIR%\runall_diffs.txt
echo ****************************** >> %TESTDIR%\runall_diffs.txt
echo. >> %TESTDIR%\runall_diffs.txt
cd LikelihoodTest
start /low /b /wait python LikelihoodTest.py
fc simulated.nex reference_output\simulated.nex >> %TESTDIR%\runall_diffs.txt
if errorlevel 1 (goto abort)
fc check.nex reference_output\check.nex >> %TESTDIR%\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

REM echo.
REM echo **********************
REM echo *** Running Phycas ***
REM echo **********************
REM echo. >> runall_diffs.txt
REM echo ********************** >> %TESTDIR%\runall_diffs.txt
REM echo *** Running Phycas *** >> %TESTDIR%\runall_diffs.txt
REM echo ********************** >> %TESTDIR%\runall_diffs.txt
REM echo. >> %TESTDIR%\runall_diffs.txt
REM cd Phycas
REM start /low /b /wait python Phycas.py
REM fc green.nex.p reference_output\green.nex.p >> %TESTDIR%\runall_diffs.txt 
REM if errorlevel 1 (goto abort)
REM fc green.nex.t reference_output\green.nex.t >> %TESTDIR%\runall_diffs.txt
REM if errorlevel 1 (goto abort)
REM cd ..

echo.
echo **************************
echo *** Running Polytomies ***
echo **************************
echo. >> runall_diffs.txt
echo ************************** >> %TESTDIR%\runall_diffs.txt
echo *** Running Polytomies *** >> %TESTDIR%\runall_diffs.txt
echo ************************** >> %TESTDIR%\runall_diffs.txt
echo. >> %TESTDIR%\runall_diffs.txt
cd Polytomies
start /low /b /wait python Polytomies.py
fc analHKY.nex.p reference_output\analHKY.nex.p >> %TESTDIR%\runall_diffs.txt
if errorlevel 1 (goto abort)
fc analHKY.nex.t reference_output\analHKY.nex.t >> %TESTDIR%\runall_diffs.txt
if errorlevel 1 (goto abort)
fc simHKY.nex reference_output\simHKY.nex >> %TESTDIR%\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

echo.
echo *************************
echo *** Running Simulator ***
echo *************************
echo. >> runall_diffs.txt
echo ************************* >> %TESTDIR%\runall_diffs.txt
echo *** Running Simulator *** >> %TESTDIR%\runall_diffs.txt
echo ************************* >> %TESTDIR%\runall_diffs.txt
echo. >> %TESTDIR%\runall_diffs.txt
cd Simulator
start /low /b /wait python Simulator.py
fc simulated.nex reference_output/simulated.nex >> %TESTDIR%\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

echo.
echo Results of tests are in the %TESTDIR%\runall_diffs.txt file
echo.
type %TESTDIR%\runall_diffs.txt
pause

goto end

:abort
echo Batch file aborted because of failed example
pause

:end
