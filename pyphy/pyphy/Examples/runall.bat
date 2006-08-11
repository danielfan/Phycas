@echo off
if exist runall_diffs.txt del runall_diffs.txt

echo ****************************
echo *** Running ExplorePrior ***
echo ****************************
echo **************************** >> runall_diffs.txt
echo *** Running ExplorePrior *** >> runall_diffs.txt
echo **************************** >> runall_diffs.txt
echo. >> runall_diffs.txt
cd ExplorePrior
start /low /b /wait python ExplorePrior.py
fc nodata.nex.p reference_output\nodata.nex.p >> ..\runall_diffs.txt 
if errorlevel 1 (goto abort)
fc nodata.nex.t reference_output\nodata.nex.t >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

echo.
echo ***************************
echo *** Running FixedParams ***
echo ***************************
echo. >> runall_diffs.txt
echo *************************** >> runall_diffs.txt
echo *** Running FixedParams *** >> runall_diffs.txt
echo *************************** >> runall_diffs.txt
echo. >> runall_diffs.txt
cd FixedParams
start /low /b /wait python FixedParams.py
fc params.p reference_output\params.p >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
fc trees.t reference_output\trees.t >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

echo.
echo ****************************
echo *** Running GelfandGhosh ***
echo ****************************
echo. >> runall_diffs.txt
echo **************************** >> runall_diffs.txt
echo *** Running GelfandGhosh *** >> runall_diffs.txt
echo **************************** >> runall_diffs.txt
echo. >> runall_diffs.txt
cd GelfandGhosh
start /low /b /wait python GelfandGhosh.py
fc ggout.txt reference_output\ggout.txt >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

echo.
echo ******************************
echo *** Running LikelihoodTest ***
echo ******************************
echo. >> runall_diffs.txt
echo ****************************** >> runall_diffs.txt
echo *** Running LikelihoodTest *** >> runall_diffs.txt
echo ****************************** >> runall_diffs.txt
echo. >> runall_diffs.txt
cd LikelihoodTest
start /low /b /wait python LikelihoodTest.py
fc simulated.nex reference_output\simulated.nex >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
fc check.nex reference_output\check.nex >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

REM echo.
REM echo **********************
REM echo *** Running Phycas ***
REM echo **********************
REM echo. >> runall_diffs.txt
REM echo ********************** >> runall_diffs.txt
REM echo *** Running Phycas *** >> runall_diffs.txt
REM echo ********************** >> runall_diffs.txt
REM echo. >> runall_diffs.txt
REM cd Phycas
REM start /low /b /wait python Phycas.py
REM fc green.nex.p reference_output\green.nex.p >> ..\runall_diffs.txt 
REM if errorlevel 1 (goto abort)
REM fc green.nex.t reference_output\green.nex.t >> ..\runall_diffs.txt
REM if errorlevel 1 (goto abort)
REM cd ..

echo.
echo **************************
echo *** Running Polytomies ***
echo **************************
echo. >> runall_diffs.txt
echo ************************** >> runall_diffs.txt
echo *** Running Polytomies *** >> runall_diffs.txt
echo ************************** >> runall_diffs.txt
echo. >> runall_diffs.txt
cd Polytomies
start /low /b /wait python Polytomies.py
fc analHKY.nex.p reference_output\analHKY.nex.p >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
fc analHKY.nex.t reference_output\analHKY.nex.t >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
fc simHKY.nex reference_output\simHKY.nex >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

echo.
echo *************************
echo *** Running Simulator ***
echo *************************
echo. >> runall_diffs.txt
echo ************************* >> runall_diffs.txt
echo *** Running Simulator *** >> runall_diffs.txt
echo ************************* >> runall_diffs.txt
echo. >> runall_diffs.txt
cd Simulator
start /low /b /wait python Simulator.py
fc simulated.nex reference_output/simulated.nex >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

goto end

:abort
echo Batch file aborted because of failed example
pause

:end
