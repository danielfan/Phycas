#! /usr/bin/env python
import os, sys
debugging = True
def mcmcOutputs(prefList):
    return [(pref + suff) for pref in prefList for suff in [".p", ".t"]]

def debug(message):
    if debugging:
        print message

def warn(message):
    print >>sys.stderr, message

def removeFilesIfTheyExist(par, fileList):
    for f in fileList:
        p = os.path.join(par, f)
        if os.path.exists(p):
            try:
                os.remove(p)
            except:
                warn(p + " not removed!")
        else:
            debug(p + " not found.")

scriptPar = os.path.split(os.path.abspath(sys.argv[0]))[0]

examplesDir = os.path.join(scriptPar, "Examples")

removeFilesIfTheyExist(scriptPar, mcmcOutputs(["runall_diffs.txt"]))
removeFilesIfTheyExist(os.path.join(examplesDir, "Data",), mcmcOutputs(["nyldna4.nex"]))
removeFilesIfTheyExist(os.path.join(examplesDir, "ExplorePrior",), mcmcOutputs(["nodata.nex"]))
removeFilesIfTheyExist(os.path.join(examplesDir, "FixedParams",), mcmcOutputs(["params"]))
ggFiles = ["ggout.txt", "simHKYg.nex"] + mcmcOutputs(["analHKY.nex", "analHKYflex.nex", "analHKYg.nex"])
removeFilesIfTheyExist(os.path.join(examplesDir, "GelfandGhosh"), ggFiles)
removeFilesIfTheyExist(os.path.join(examplesDir, "LikelihoodTest",), ["simulated.nex", "check.nex"])
removeFilesIfTheyExist(os.path.join(examplesDir, "Polytomies",), ["simHKY.nex"] + mcmcOutputs(["analHKY.nex"]))
removeFilesIfTheyExist(os.path.join(examplesDir, "Simulator",), ["simulated.nex"])

