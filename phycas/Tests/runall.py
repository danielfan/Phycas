import os, sys, re
import difflib
debugging = False

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

def writeHeader(out, name):
    s = "*** Running " + name + " ***"
    border = "*" * len(s)
    out.write("%s\n%s\n%s\n \n" % (border, s, border))

_diff_key_pat = re.compile(r"[-*]{3}\s+(\d+)(,?\d*)\s+[-*]{3}")
_missig_key_pat = re.compile(r"[-]{3}\s+(\d+)(,?\d*)\s+[-]{3}")

def parse_diffs(diff_lines):
    d = {}
    k, t, curr = None, [], []
    for line in diff_lines:
        if line.strip() == "***************":
            if curr:
                t.append(curr)
            if k is not None:
                d[k] = tuple(t)
            k, t, curr = None, [], []
        else:
            m = _diff_key_pat.match(line)
            if m:
                if k is None:
                    k = int(m.group(1))
                else:
                    t.append(curr)
                t.append([line])
                curr = []
            else:
                curr.append(line)
    if k is not None:
        if curr:
            t.append(curr)
        d[k] = tuple(t)
    return d

def runDiff(a, b, outFile, acceptable_diffs=()):
    fa = open(a, "rU")
    fb = open(b, "rU")
    all_diffs = list(difflib.context_diff(fa.readlines(), fb.readlines(), n=0))
    print "\n".join(all_diffs)
    if all_diffs:
        all_diffs = parse_diffs(all_diffs)
        for k in acceptable_diffs:
            try:
                del all_diffs[k]
            except: 
                pass
    if all_diffs:
        l = []
        for k, vi in all_diffs.iteritems():
            t = []
            for v in vi:
                t.extend(v)
            l.append((k, t))
        l.sort()
        nd = []
        for k, lines in l:
            nd.extend(lines)
        s = "".join(nd)
        fout = open(outFile, "a")
        fout.write("Diffs in file %s:\n" % a)
        fout.write(s)
        if len(s) < 5000:
            t = ":\n%s\n(diffs also stored at the end of %s)." % (s, outFile)
        else:
            t = ". See the end of the %s." % (outFile)
        fout.close()
        sys.exit("Script aborted because of failed example.\nOutput differed from the expected output%s\n" % t)
    fa.close()
    fb.close()
    
def runTest(outFile, name, results):
    print("\n")
    print("________________________________________")
    print("|\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/|")
    print("| Running test in " + name)
    print("|\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/|")
    print("----------------------------------------")
    print("\n")
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
        acceptable_diff_f = os.path.join('acceptable_diff', f)
        if os.path.exists(acceptable_diff_f):
            # right now we simpli store a line numbers in "acceptable_diff_f/filename"
            #   for lines that can generate a diff
            acceptable_diff = [int(i.strip()) for i in open(acceptable_diff_f, "rU")]
        else:
            acceptable_diff = ()
        runDiff(f, os.path.join('reference_output', f), outFile, acceptable_diff)
    os.chdir(prevDir)

if __name__ == '__main__':
    if 'python_d' in sys.executable:
        debugging = True
    scriptPar = os.path.split(os.path.abspath(sys.argv[0]))[0]
    outFile = os.path.join(scriptPar, "diffs.txt")
    removeFilesIfTheyExist(scriptPar, [outFile])

    # os.chdir(os.path.join(os.path.split(scriptPar)[0], "Examples"))
    os.chdir(scriptPar)
    runTest(outFile, "Sumt", ["trees.tre","splits.pdf","logfile.txt"])
    runTest(outFile, "PathSampling", ["params.p", "trees.t"])
    runTest(outFile, "LikelihoodTest", ["simulated.nex", "check.nex"])
    #runTest(outFile, "ExplorePrior", mcmcOutputs(["nodata.nex"]))
    runTest(outFile, "FixedParams", ["fixed.p", "fixed.t"])
    #runTest(outFile, "FixedTopology", ["fixdtree.p", "fixdtree.t", "simulated.nex"])
    runTest(outFile, "GTRTest", ["gtr_test.p", "gtr_test.t"])
    runTest(outFile, "Simulator", ["simulated.nex"])
    runTest(outFile, "SplitTest", ["out.txt"])
    runTest(outFile, "PDFTree", ["test.pdf"])
    # note: should add trees.pdf to list for SumT, but slight rounding differences
    # cause PDF files to be different, and haven't been able to figure out
    # why the rounding should be different between Intel Mac and Intel PC!

    # Still need to get these working...
    #runTest(outFile, "SAMCTest", ["cf.txt", "samc_output.p", "samc_output.t"])
    #runTest(outFile, "GelfandGhosh", ["ggout.txt", "analHKY.nex.p", "analHKY.nex.t", "analHKYflex.nex.p", "analHKYflex.nex.t", "analHKYg.nex.p", "analHKYg.nex.t"])
    #runTest(outFile, "Polytomies", ["simHKY.nex"] + mcmcOutputs(["HKYpolytomy"]))

    if sys.platform == 'win32':
        raw_input('Press any key to quit')
