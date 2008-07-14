import os, sys, subprocess

_test_dir = None

def getPhycasTestDir():
    global _test_dir
    if _test_dir is None:
        d = os.path.dirname(os.path.dirname(__file__))
        _test_dir = os.path.join(d, "Tests")
    return _test_dir

def getPhycasTestDataDir():
    return os.path.join(getPhycasTestDir(),"Data")

def getPhycasTestData(filen):
    return os.path.join(getPhycasTestDataDir(), filen)

def _runRegressionTests(out):
    d = getPhycasTestDir()
    r = os.path.join(d, "runall.py")
    spawnPython(r)
    r = os.path.join(d, "doctestall.py")
    spawnPython(r)
    out.info("\nAll tests passed.")
    
        
def spawnPython(f):
    if not os.path.exists(f):
        raise RuntimeError('The python script "%s" does not exist' % f)
    retcode = subprocess.call([sys.executable, f])
    if retcode < 0:
        raise RuntimeError("python execution of %s was terminated by signal %s" % (f, str(-retcode)))
    elif retcode > 0:
        raise RuntimeError("python execution of %s was failed with code %s" % (f, str(retcode)))
