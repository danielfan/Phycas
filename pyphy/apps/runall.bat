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
python ExplorePrior.py
fc nodata.nex.p refExplorePrior\nodata.nex.p >> ..\runall_diffs.txt 
if errorlevel 1 (goto abort)
fc nodata.nex.t refExplorePrior\nodata.nex.t >> ..\runall_diffs.txt
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
python FixedParams.py
fc params.p ref_output\params.p >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
fc trees.t ref_output\trees.t >> ..\runall_diffs.txt
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
python GelfandGhosh.py
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
python LikelihoodTest.py
fc simulated.nex refLikelihoodTest\simulated.nex >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
fc check.nex refLikelihoodTest\check.nex >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

echo.
echo **************************
echo *** Running MCMCSimple ***
echo **************************
echo. >> runall_diffs.txt
echo ************************** >> runall_diffs.txt
echo *** Running MCMCSimple *** >> runall_diffs.txt
echo ************************** >> runall_diffs.txt
echo. >> runall_diffs.txt
cd MCMCSimple
python MCMCSimple.py
fc nyldna4.nex.p refHKYhyper\nyldna4.nex.p >> ..\runall_diffs.txt 
if errorlevel 1 (goto abort)
fc nyldna4.nex.t refHKYhyper\nyldna4.nex.t >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

echo.
echo **********************
echo *** Running Phycas ***
echo **********************
echo. >> runall_diffs.txt
echo ********************** >> runall_diffs.txt
echo *** Running Phycas *** >> runall_diffs.txt
echo ********************** >> runall_diffs.txt
echo. >> runall_diffs.txt
cd Phycas
python Phycas.py
fc nyldna4.nex.p refHKYhyper\nyldna4.nex.p >> ..\runall_diffs.txt 
if errorlevel 1 (goto abort)
fc nyldna4.nex.t refHKYhyper\nyldna4.nex.t >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

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
python Polytomies.py
fc analHKY.nex.p refHKY\analHKY.nex.p >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
fc analHKY.nex.t refHKY\analHKY.nex.t >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
fc simHKY.nex refHKY\simHKY.nex >> ..\runall_diffs.txt
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
python Simulator.py
fc simulated.nex refHKY/simulated.nex >> ..\runall_diffs.txt
if errorlevel 1 (goto abort)
cd ..

goto end

:abort
echo Batch file aborted because of failed example
pause

:end
