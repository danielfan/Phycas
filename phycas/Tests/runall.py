import os, sys

debugging = False
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
    interpreter = 'python'
    if sys.platform == 'win32' and debugging:
        interpreter = 'python_d'
    if os.system('%s %s.py' % (interpreter,name)) != 0:
        sys.exit("Script aborted because of failed example")
    for f in results:
        runDiff(f, os.path.join('reference_output', f), outFile)
    os.chdir(prevDir)
    
scriptPar = os.path.split(os.path.abspath(sys.argv[0]))[0]
outFile = os.path.join(scriptPar, "runall_diffs.txt")
removeFilesIfTheyExist(scriptPar, [outFile])

# os.chdir(os.path.join(os.path.split(scriptPar)[0], "Examples"))

runTest(outFile, "ExplorePrior", mcmcOutputs(["nodata.nex"]))
runTest(outFile, "FixedParams", ["fixed.p", "fixed.t"])
#runTest(outFile, "GelfandGhosh", ["ggout.txt", "analHKY.nex.p", "analHKY.nex.t", "analHKYflex.nex.p", "analHKYflex.nex.t", "analHKYg.nex.p", "analHKYg.nex.t"])
runTest(outFile, "GTRTest", ["gtr_test.p", "gtr_test.t"])
#runTest(outFile, "LikelihoodTest", ["simulated.nex", "check.nex"])
#runTest(outFile, "Polytomies", ["simHKY.nex"] + mcmcOutputs(["HKYpolytomy"]))
runTest(outFile, "Simulator", ["simulated.nex"])

if sys.platform == 'win32':
    raw_input('Press any key to quit')
