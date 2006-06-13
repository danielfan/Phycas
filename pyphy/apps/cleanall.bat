@echo off
if exist runall_diffs.txt del /Q runall_diffs.txt

if exist ExplorePrior\nodata.nex.p del /Q ExplorePrior\nodata.nex.p
if exist ExplorePrior\nodata.nex.t del /Q ExplorePrior\nodata.nex.t

if exist FixedParams\params.p del /Q FixedParams\params.p
if exist FixedParams\trees.t del /Q FixedParams\trees.t

if exist GelfandGhosh\ggout.txt del /Q GelfandGhosh\ggout.txt
if exist GelfandGhosh\analHKY.nex.p del /Q GelfandGhosh\analHKY.nex.p
if exist GelfandGhosh\analHKY.nex.t del /Q GelfandGhosh\analHKY.nex.t
if exist GelfandGhosh\analHKYflex.nex.p del /Q GelfandGhosh\analHKYflex.nex.p
if exist GelfandGhosh\analHKYflex.nex.t del /Q GelfandGhosh\analHKYflex.nex.t
if exist GelfandGhosh\analHKYg.nex.p del /Q GelfandGhosh\analHKYg.nex.p
if exist GelfandGhosh\analHKYg.nex.t del /Q GelfandGhosh\analHKYg.nex.t
if exist GelfandGhosh\simHKYg.nex del /Q GelfandGhosh\simHKYg.nex

if exist LikelihoodTest\simulated.nex del /Q LikelihoodTest\simulated.nex
if exist LikelihoodTest\check.nex del /Q LikelihoodTest\check.nex
                       
if exist MCMCSimple\nyldna4.nex.p del /Q MCMCSimple\nyldna4.nex.p
if exist MCMCSimple\nyldna4.nex.t del /Q MCMCSimple\nyldna4.nex.t

if exist Phycas\green.nex.p del /Q Phycas\green.nex.p
if exist Phycas\green.nex.t del /Q Phycas\green.nex.t

if exist Polytomies\analHKY.nex.p del /Q Polytomies\analHKY.nex.p
if exist Polytomies\analHKY.nex.t del /Q Polytomies\analHKY.nex.t
if exist Polytomies\simHKY.nex del /Q Polytomies\simHKY.nex

if exist Simulator\simulated.nex del /Q Simulator\simulated.nex
