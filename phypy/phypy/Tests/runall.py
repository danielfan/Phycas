import os, sys

debugging = True
diffUtility = sys.platform == 'win32' and 'fc' or 'diff'

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
            finally:
                warn(p + " not removed!")
        else:
            debug(p + " not found.")

def writeHeader(out, name):
    s = "*** Running " + name + " ***"
    border = "*" * len(s)
    out.write("%s\n%s\n%s\n \n" % (border, s, border))

def runDiff(a, b, outFile):
    if os.system('%s %s %s >> %s' % (diffUtility, a, b, outFile)) != 0:
         sys.exit("Script aborted because of failed example.\nOutput differed from the expected output.  See the differences at the end of %s" % outFile)
    
def runTest(outFile, name, results):
    outStream = open(outFile, 'a')
    writeHeader(outStream, name)
    outStream.close()
    prevDir = os.path.abspath(os.curdir)
    os.chdir(name)
    if os.system('python %s.py' % name) != 0:
        sys.exit("Script aborted because of failed example")
    for f in results:
        runDiff(f, os.path.join('reference_output', f), outFile)
    os.chdir(prevDir)
    
scriptPar = os.path.split(os.path.abspath(sys.argv[0]))[0]
outFile = os.path.join(scriptPar, "runall_diffs.txt")
removeFilesIfTheyExist(scriptPar, [outFile])

# os.chdir(os.path.join(os.path.split(scriptPar)[0], "Examples"))

runTest(outFile, "ExplorePrior", mcmcOutputs(["nodata.nex"]))
runTest(outFile, "FixedParams", ["params.p", "trees.t"])
runTest(outFile, "GelfandGhosh", ["ggout.txt"])
runTest(outFile, "LikelihoodTest", ["simulated.nex", "check.nex"])
runTest(outFile, "Polytomies", ["simHKY.nex"] + mcmcOutputs(["analHKY.nex"]))
runTest(outFile, "Simulator", ["simulated.nex"])
#runTest(outFile, "Phycas", mcmcOutputs(["green.nex"]))

raw_input('Press any key to quit')
    