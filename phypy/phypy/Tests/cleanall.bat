@echo off
if exist runall_diffs.txt del /Q runall_diffs.txt

if exist ..\Examples\Data\nyldna4.nex.p del /Q ..\Examples\Data\nyldna4.nex.p
if exist ..\Examples\Data\nyldna4.nex.t del /Q ..\Examples\Data\nyldna4.nex.t

if exist ..\Examples\ExplorePrior\nodata.nex.p del /Q ..\Examples\ExplorePrior\nodata.nex.p
if exist ..\Examples\ExplorePrior\nodata.nex.t del /Q ..\Examples\ExplorePrior\nodata.nex.t

if exist ..\Examples\FixedParams\params.p del /Q ..\Examples\FixedParams\params.p
if exist ..\Examples\FixedParams\trees.t del /Q ..\Examples\FixedParams\trees.t

if exist ..\Examples\GelfandGhosh\ggout.txt del /Q ..\Examples\GelfandGhosh\ggout.txt
if exist ..\Examples\GelfandGhosh\analHKY.nex.p del /Q ..\Examples\GelfandGhosh\analHKY.nex.p
if exist ..\Examples\GelfandGhosh\analHKY.nex.t del /Q ..\Examples\GelfandGhosh\analHKY.nex.t
if exist ..\Examples\GelfandGhosh\analHKYflex.nex.p del /Q ..\Examples\GelfandGhosh\analHKYflex.nex.p
if exist ..\Examples\GelfandGhosh\analHKYflex.nex.t del /Q ..\Examples\GelfandGhosh\analHKYflex.nex.t
if exist ..\Examples\GelfandGhosh\analHKYg.nex.p del /Q ..\Examples\GelfandGhosh\analHKYg.nex.p
if exist ..\Examples\GelfandGhosh\analHKYg.nex.t del /Q ..\Examples\GelfandGhosh\analHKYg.nex.t
if exist ..\Examples\GelfandGhosh\simHKYg.nex del /Q ..\Examples\GelfandGhosh\simHKYg.nex

if exist ..\Examples\LikelihoodTest\simulated.nex del /Q ..\Examples\LikelihoodTest\simulated.nex
if exist ..\Examples\LikelihoodTest\check.nex del /Q ..\Examples\LikelihoodTest\check.nex
                       
if exist ..\Examples\Polytomies\analHKY.nex.p del /Q ..\Examples\Polytomies\analHKY.nex.p
if exist ..\Examples\Polytomies\analHKY.nex.t del /Q ..\Examples\Polytomies\analHKY.nex.t
if exist ..\Examples\Polytomies\simHKY.nex del /Q ..\Examples\Polytomies\simHKY.nex

if exist ..\Examples\Simulator\simulated.nex del /Q ..\Examples\Simulator\simulated.nex
